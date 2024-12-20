# Get a lightcurve from a folder of FITs files

import sep_pjw as sep
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time
from os import walk
import math
import numpy as np
import random
import sys
import time
random.seed(int(time.time()))

# Whether or not to use command line args
USE_CLI_ARGS = True
# other general config
FILE_FOLDER_PATH = "../science-data/corot-1b/"
STAR_NAME = "CoRoT-1 b"
FILTER = "V"
MAG_MODE = True # whether or not to use magnitudes instead of raw flux ratio
RANDOM_INSPECTION = False # randomly inspects a file so you can manually check the apertures
VERBOSE = True # extra logging info in console
EXPORT_CSV = True # whether or not to export photometry to csv
EXPORT_PATH = f"./{STAR_NAME.lower().replace(" ", "-")}-{FILTER.lower()}.csv"

# ~~~ Choose a target ~~~ #

# BW Cru
# TARGET_RA = "12 53 57.5386794336"
# TARGET_DEC = "-60 24 58.092625152"

# CV Cru
# TARGET_RA = "12 53 46.9935752640"
# TARGET_DEC = "-60 18 35.656553808"

TARGET_RA = ""
TARGET_DEC = ""

# if no WCS, insert the initial x, y pixel values of the target star here
TARGET_PIX = [0, 0]

# ~~~ Choose a comp(s) ~~~ #

# HD 312084
# COMPS_RA = ["12 54 14.556"]
# COMPS_DEC = ["-60 21 47.28"]
# COMPS_MAG_B = [10.92]
# COMPS_MAG_V = [10.79]

# HD 312080
# COMPS_RA = ["12 53 14.258"]
# COMPS_DEC = ["-60 27 38.129"]
# COMPS_MAG_B = [10.23]
# COMPS_MAG_V = [10.10]

# HD 111990
# COMPS_RA = ["12 53 59.7993"]
# COMPS_DEC = ["-60 20 07.5103"]
# COMPS_MAG_B = [6.95]
# COMPS_MAG_V = [6.78]

# HD 111714, HD 111699, HD 111916
# COMPS_RA = ["12 52 01.704", "12 51 52.435", "12 53 26.557"]
# COMPS_DEC = ["-60 33 01.44", "-60 24 21.83", "-60 12 00.37"]
# COMPS_MAG_B = [9.465, 9.86, 9.549]
# COMPS_MAG_V = [9.339, 9.754, 9.313]
# COMPS_MAG_up = [10.46, 10.86, 10.57]

COMPS_RA = [""]
COMPS_DEC = [""]
COMPS_MAG = [0]

# if no WCS, insert the initial x, y pixel values of the comp star(s) here
COMPS_PIX = [[0, 0]]


if USE_CLI_ARGS:
	TARGET_RA, TARGET_DEC = "", ""
	TARGET_PIX = [-1, -1]
	COMPS_RA, COMPS_DEC = [""], [""]
	COMPS_PIX = [[-1, -1]]
	VERBOSE, MAG_MODE, RANDOM_INSPECTION, EXPORT_CSV = (False,)*4
	req_var_set = 7 # required variable sets
	i = 1
	while i < len(sys.argv):
		req_var_set -= 1
		t_f = sys.argv[i] # this_flag
		if t_f == "--path": # data folder path
			FILE_FOLDER_PATH = sys.argv[i+1]
		elif t_f == "--filter": # filter of the obs
			FILTER = sys.argv[i+1]
		elif t_f == "--tname": # target star name
			STAR_NAME = sys.argv[i+1]
		elif t_f == "--tra": # target star RA
			TARGET_RA = sys.argv[i+1]
		elif t_f == "--tdec": # target star Dec
			TARGET_DEC = sys.argv[i+1]
		elif t_f == "--tpix": # target pixel values, separated by space
			TARGET_PIX = map(float, sys.argv[i+1].split(","))
		elif t_f == "--cra": # comp star RA
			COMPS_RA = [sys.argv[i+1]]
		elif t_f == "--cdec": # comp star Dec
			COMPS_DEC = [sys.argv[i+1]]
		elif t_f == "--cpix": # comparison star pixel values, separated by space
			COMPS_PIX = [map(float, sys.argv[i+1].split(","))]
		elif t_f == "-i": # whether or not to do inspection
			RANDOM_INSPECTION = True
			req_var_set += 1
			i -= 1
		elif t_f == "--cmag": # comp star mag
			MAG_MODE = True
			COMPS_MAG = [float(sys.argv[i+1])]
			req_var_set += 1
		elif t_f == "-v": # verbose
			VERBOSE = True
			req_var_set += 1
		elif t_f == "-e": # export data to CSV
			EXPORT_CSV = True
			req_var_set += 1
		i += 2
	if req_var_set > 0:
		print("Not all required flags have been set!")
		exit(1)
	elif TARGET_RA == "" and TARGET_DEC == "" and TARGET_PIX == [-1, -1]:
		print("Target star position has not been set")
		exit(1)
	elif COMPS_RA == [""] and COMPS_DEC == [""] and COMPS_PIX == [[-1, -1]]:
		print("Comparison star position has not been set")
		exit(1)

print("Initializing SEP photometer with these settings:")
print("File folder path:", FILE_FOLDER_PATH)
print("Filter:", FILTER)
print("Target star RA:", TARGET_RA)
print("Target star Dec:", TARGET_DEC)
print("Comp star(s) RA:", COMPS_RA)
print("Comp star(s) Dec:", COMPS_DEC)
print("Inspection:", RANDOM_INSPECTION)
print("Magnitude mode:", MAG_MODE)
if MAG_MODE == True:
	print("Comp star(s) Magnitudes", COMPS_MAG)
print()

# gets the average background in an annulus with inner radius r and outer radius r + dr
def sky_bg_from_annulus(data, cx, cy, r, dr):
	small_mask = np.zeros(data.shape, dtype=bool)
	big_mask = np.zeros(data.shape, dtype=bool)
	sep.mask_ellipse(small_mask, [cx], [cy], [r], [r], [0])
	sep.mask_ellipse(big_mask, [cx], [cy], [r+dr], [r+dr], [0])
	annulus_mask = np.logical_xor(small_mask, big_mask)
	annulus_pixels = []
	for y in range(int(cy - (r + dr) - 1), int(cy + (r + dr) + 2)):
		for x in range(int(cx - (r + dr) - 1), int(cx + (r + dr) + 2)):
			if annulus_mask[y][x]:
				annulus_pixels.append(data[y][x])
	cutoff = np.percentile(annulus_pixels, 99)
	for i in range(len(annulus_pixels)-1, -1, -1):
		if annulus_pixels[i] > cutoff:
			del annulus_pixels[i]
	return np.mean(annulus_pixels)


TARGET_SKYCOORD = SkyCoord(ra=TARGET_RA, dec=TARGET_DEC, unit=(u.hourangle, u.deg))
COMPS_SKYCOORD = []
for i in range(len(COMPS_RA)):
	COMPS_SKYCOORD.append(SkyCoord(ra=COMPS_RA[i], dec=COMPS_DEC[i], unit=(u.hourangle, u.deg)))

# arrays for exporting later
MJD = []
FLUX_RATIO = []
FLUX_RATIO_ERR = []

fits_filenames = []
for (_, _, filenames) in walk(FILE_FOLDER_PATH):
    fits_filenames.extend(filenames)
    break

# if RANDOM_INSPECTION == True:
# 	RANDOM_INSPECTION = FILE_FOLDER_PATH + random.choice(fits_filenames)

START_TIME = int(time.time())
RAW_TARGET_FLUX = []
RAW_COMPS_FLUX = []
for file_name in fits_filenames:
	file_path = FILE_FOLDER_PATH + file_name
	with fits.open(file_path) as hdul:
		real_hdul = hdul[1] if file_name.endswith(".fits.fz") else hdul[0]
		if not "filter" in real_hdul.header:
			continue
		if real_hdul.header["filter"].strip() != FILTER or "-e00." in file_name:
			continue
		print("Opening file path", file_path)
		use_wcs = True
		this_wcs = WCS(real_hdul.header)
		if this_wcs.wcs.ctype[0] != "":
			target_x, target_y = this_wcs.world_to_pixel(TARGET_SKYCOORD)
		else:
			use_wcs = False
			target_x, target_y = TARGET_PIX
		fits_data = real_hdul.data.byteswap().newbyteorder()
		if file_name.endswith(".fits.fz"):
			fits_data = real_hdul.data
		sep_bkg = sep.Background(fits_data)
		sep_bkg_err = sep_bkg.globalrms

		target_x, target_y, _ = sep.winpos(fits_data, [target_x], [target_y], [3.33])
		target_x = target_x[0]
		target_y = target_y[0]

		# optimize radius from flux_radius and kron_radius
		loc_sky_bg = sky_bg_from_annulus(fits_data, target_x, target_y, 12, 8)
		init_rad = sep.flux_radius(fits_data - loc_sky_bg, [target_x], [target_y], [12], 0.66)[0][0]
		kron_rad, kr_flag = sep.kron_radius(fits_data, [target_x], [target_y], [2*init_rad], [2*init_rad], [0], [6])

		if kr_flag[0] != 0:
			print("\nSEP raised a non-zero flag for target star kron radius: ", kr_flag[0])
			exit(kr_flag[0])
		kron_rad = kron_rad[0] * 1.35

		target_flux, target_err, t_flag = sep.sum_circle(fits_data, [target_x], [target_y], [kron_rad],
			err=sep_bkg_err, bkgann=(kron_rad+4, kron_rad+12))
		if t_flag[0] != 0:
			print("\nSEP raised a non-zero flag for targetg star photometry:", t_flag[0])
			exit(t_flag[0])
		target_flux = target_flux[0]
		target_err = target_err[0]
		if VERBOSE == True:
			print("Target X Y coords and radius", target_x, target_y, kron_rad)

		target_mag = 0 # only for mag mode
		total_comp_err = 0
		total_comp_flux = 0 # total comp flux is used for both ADU counts and magnitude
		mat_patches = [] # matplotlib patches
		for i in range(max(len(COMPS_SKYCOORD),len(COMPS_PIX))):
			if use_wcs:
				comp_coord = COMPS_SKYCOORD[i]
				comp_x, comp_y = this_wcs.world_to_pixel(comp_coord)
			else:
				comp_x, comp_y = COMPS_PIX[i]
			comp_x, comp_y, _ = sep.winpos(fits_data, [comp_x], [comp_y], [3.33])
			comp_x = comp_x[0]
			comp_y = comp_y[0]

			loc_sky_bg = sky_bg_from_annulus(fits_data, comp_x, comp_y, 12, 8)
			init_rad = sep.flux_radius(fits_data - loc_sky_bg, [comp_x], [comp_y], [12], 0.66)[0][0]
			comp_krad, ckr_flag = sep.kron_radius(fits_data, [comp_x], [comp_y], [2*init_rad], [2*init_rad], [0], [6])

			if ckr_flag[0] != 0:
				print(f"\nSEP raised a non-zero flag for comp star {i} kron radius: ", ckr_flag[0])
				exit(ckr_flag[0])
			comp_krad = comp_krad[0] * 1.35

			comp_flux, comp_err, c_flag = sep.sum_circle(fits_data, [comp_x], [comp_y], [comp_krad],
				err=sep_bkg_err, bkgann=(comp_krad+4, comp_krad+12))
			if c_flag[0] !=0:
				print(f"\nSEP raised a non-zero flag with comp star {i+1} aperture photometry:", c_flag[0])
				exit(t_flag[0])
			comp_flux = comp_flux[0]
			comp_err = comp_err[0]

			if VERBOSE == True:
				print("Comp X Y coords and radius", comp_x, comp_y, comp_krad)

			circ2 = Circle((comp_x, comp_y), comp_krad, color="r", fill=False)
			mat_patches.append(circ2)
			if MAG_MODE == True:
				this_comp_mag = COMPS_MAG[i]
				this_target_mag = this_comp_mag - 2.5 * math.log10(target_flux / comp_flux)
				target_mag += this_target_mag / max(len(COMPS_SKYCOORD), len(COMPS_PIX))
			total_comp_flux += comp_flux
			total_comp_err += comp_err
			COMPS_PIX[i] = [comp_x, comp_y]
		total_target_err = math.sqrt((target_err / target_flux) ** 2 + (total_comp_err / total_comp_flux)**2)
		ratio = target_flux / total_comp_flux
		if MAG_MODE == True:
			ratio = target_mag
			total_target_err = 2.5 / math.log(10) * total_target_err

		if RANDOM_INSPECTION == True and random.uniform(0, 1) < 0.02:
			fig = plt.figure(figsize=(8, 6))
			ax = fig.add_subplot()
			ax.set_title(f"Sample image of target star (blue) and comps (red)")
			im1 = ax.imshow(fits_data, interpolation="nearest", vmin=0.0, vmax=4000.0, cmap="gray", origin="lower")
			fig.colorbar(im1)
			circ = Circle((target_x, target_y), kron_rad, color="b", fill=False)
			ax.add_patch(circ)
			for p in mat_patches:
				ax.add_patch(p)
			plt.show()
			RANDOM_INSPECTION = False

		mid_photo_time = 0
		if "MJD-OBS" in real_hdul.header:
			mid_photo_time = real_hdul.header["MJD-OBS"]
		elif "DATE-AVG" in real_hdul.header:
			mid_photo_time = Time(real_hdul.header["DATE-AVG"]).jd
		elif "DATE-OBS" in real_hdul.header:
			mid_photo_time = Time(real_hdul.header["DATE-OBS"]).jd
		else:
			print("\nNo suitable date header found!!! :(")
			exit()
		MJD.append(mid_photo_time)
		FLUX_RATIO.append(ratio)
		FLUX_RATIO_ERR.append(total_target_err)
		TARGET_PIX = [target_x, target_y]
		RAW_TARGET_FLUX.append(target_flux)
		RAW_COMPS_FLUX.append(total_comp_flux)
		print()

# sort all the data points by MJD so it can be displayed easily
COMBINED_DATA = []
MIN_MJD = min(MJD)
for i in range(len(MJD)):
	COMBINED_DATA.append([MJD[i] - MIN_MJD, FLUX_RATIO[i], FLUX_RATIO_ERR[i]])

MJD = []
FLUX_RATIO = []
FLUX_RATIO_ERR = []
COMBINED_DATA = sorted(COMBINED_DATA, key=(lambda x: x[0]))
if EXPORT_CSV == True:
	with open(EXPORT_PATH, "a") as data_file:
		data_title = "Magnitude" if MAG_MODE == True else "Flux ratio"
		data_file.write(f"MJD,{data_title},Uncertainty\n")
for data_set in COMBINED_DATA:
	MJD.append(data_set[0])
	FLUX_RATIO.append(data_set[1])
	FLUX_RATIO_ERR.append(data_set[2])
	if EXPORT_CSV == True:
		with open(EXPORT_PATH, "a") as data_file:
			data_file.write("{},{},{}\n".format(data_set[0], data_set[1], data_set[2]))

the_s = "s" if len(COMPS_SKYCOORD) > 1 else ""
ylabel = f"Flux ratio (target / comp{the_s})" if MAG_MODE == False else "Mangitude"

fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot()
ax.set_title(f"Time-series flux of {STAR_NAME} in the filter {FILTER}")
ax.errorbar(MJD, FLUX_RATIO, yerr=FLUX_RATIO_ERR, linewidth=0, elinewidth=1.5, ecolor="#333", alpha=0.5)
ax.scatter(MJD, FLUX_RATIO, marker="o", color="green", s=15)
ax.set_xlabel("Time (Days)")
ax.set_ylabel(ylabel)

END_TIME = int(time.time())
TOTAL_TIME = END_TIME - START_TIME
print(f"Complete! That took {TOTAL_TIME} seconds")

plt.show()

# debug plots showing raw flux of target and comp stars
# fig, axs = plt.subplots(2)
# fig.suptitle("Raw-Flux vs Time of Target Star, Comp Star")
# axs[0].scatter(MJD, RAW_TARGET_FLUX, marker="o", color="purple", s=15)
# axs[1].scatter(MJD, RAW_COMPS_FLUX, marker="o", color="blue", s=15)
# plt.show()
