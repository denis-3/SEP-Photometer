# Get a lightcurve from a folder of FITs files

import sep
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

# WASP 164 b
# TARGET_RA = "22 59 29.67"
# TARGET_DEC = "-60 26 52.26"

# WASP 164 b Long period variable - ASAS J225857-6028.8
# TARGET_RA = "22 58 57.334"
# TARGET_DEC = "-60 28 48.72"

# CoRoT-1 b
# TARGET_RA = "06 48 19.172"
# TARGET_DEC = "-03 06 07.86"

# V0982 Mon (near CoRoT-1 b)
TARGET_RA = "06 48 47.168"
TARGET_DEC = "-02 53 53.26"

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

# For WASP-164 b
# COMPS_RA = ["22 59 48.334"]#, "22 59 54.413", "23 00 36"]
# COMPS_DEC = ["-60 30 59.52"]#, "-60 31 18.61", "-60 21 58.18"]
# COMPS_MAG = [11.767]

# For CoRoT-1 b
COMPS_RA = ["06 48 22.111"]
COMPS_DEC = ["-03 05 15.78"]
COMPS_MAG = [13.4]


if USE_CLI_ARGS:
	VERBOSE = False
	MAG_MODE = False
	RANDOM_INSPECTION = False
	EXPORT_CSV = False
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
		elif t_f == "--cra": # comp star RA
			COMPS_RA = [sys.argv[i+1]]
		elif t_f == "--cdec": # comp star Dec
			COMPS_DEC = [sys.argv[i+1]]
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
	cutoff = np.percentile(annulus_pixels, 95)
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

for file_name in fits_filenames:
	file_path = FILE_FOLDER_PATH + file_name
	with fits.open(file_path) as hdul:
		real_hdul = hdul[1] if file_name.endswith(".fits.fz") else hdul[0]
		if not "filter" in real_hdul.header:
			continue
		if real_hdul.header["filter"].strip() != FILTER or "-e00." in file_name:
			continue
		print("Opening file path", file_path)
		this_wcs = WCS(real_hdul.header)
		target_x, target_y = this_wcs.world_to_pixel(TARGET_SKYCOORD)
		fits_data = real_hdul.data.byteswap().newbyteorder()
		if file_name.endswith(".fits.fz"):
			fits_data = real_hdul.data
		# bkg_sep = sep.Background(fits_data)
		# data_minus_bg = fits_data - bkg_sep
		data_minus_bg = fits_data

		target_x, target_y, _ = sep.winpos(data_minus_bg, [target_x], [target_y], [3.33])
		target_x = target_x[0]
		target_y = target_y[0]

		# local target background
		loc_tgt_bg = sky_bg_from_annulus(data_minus_bg, target_x, target_y, 14, 10)
		loc_data_minus_bg = data_minus_bg - loc_tgt_bg
		kron_rad = sep.flux_radius(loc_data_minus_bg, [target_x], [target_y], [11], 0.9)[0][0]

		target_flux, fluxerr, _ = sep.sum_circle(loc_data_minus_bg, [target_x], [target_y], [kron_rad], subpix=10)
		if VERBOSE == True:
			print("Target X Y coords and radius", target_x, target_y, kron_rad)

		# kron_rad = sep.kron_radius(data_minus_bg, [target_x], [target_y], 20, 20, 0, 6.0)[0][0] * 1.4
		# target_flux, fluxerr, _ = sep.sum_circle(data_minus_bg, [target_x], [target_y], kron_rad[0], subpix=1)

		total_comp_err = 0
		total_comp_flux = 0 # total comp flux is used for both ADU counts and magnitude
		target_mag = 0 # only for mag mode
		mat_patches = [] # matplotlib patches
		for i in range(len(COMPS_SKYCOORD)):
			comp_coord = COMPS_SKYCOORD[i]
			comp_x, comp_y = this_wcs.world_to_pixel(comp_coord)
			comp_x, comp_y, _ = sep.winpos(data_minus_bg, [comp_x], [comp_y], [3.33])
			comp_x = comp_x[0]
			comp_y = comp_y[0]

			loc_cmp_bg = sky_bg_from_annulus(data_minus_bg, comp_x, comp_y, 14, 10)
			loc_data_minus_bg = data_minus_bg - loc_cmp_bg

			comp_krad = sep.flux_radius(loc_data_minus_bg, [comp_x], [comp_y], [11], 0.9)[0][0]
			comp_flux, comp_ferr, _ = sep.sum_circle(loc_data_minus_bg, [comp_x], [comp_y], [comp_krad], subpix=10)
			# comp_krad = sep.kron_radius(data_minus_bg, [comp_x], [comp_y], 20, 20, 0, 6.0)[0][0] * 1.4
			# comp_flux, comp_ferr, _ = sep.sum_circle(data_minus_bg, [comp_x], [comp_y], comp_krad, subpix=1)
			if VERBOSE == True:
				print("Comp X Y coords and radius", comp_x, comp_y, comp_krad)

			circ2 = Circle((comp_x, comp_y), comp_krad, color="r", fill=False)
			mat_patches.append(circ2)
			if MAG_MODE == True:
				this_comp_mag = COMPS_MAG[i]
				this_target_mag = this_comp_mag - 2.5 * math.log10(target_flux[0] / comp_flux[0])
				target_mag += this_target_mag / len(COMPS_RA)
				total_comp_flux += this_comp_mag / len(COMPS_RA)
				total_comp_err = 0
			else:
				total_comp_flux += comp_flux[0]
				total_comp_err += comp_ferr[0]
		if MAG_MODE == True:
			ratio = target_mag
		else:
			ratio = target_flux[0] / total_comp_flux
		# total_err = (target_flux[0] + fluxerr[0]) / (total_comp_flux + total_comp_err)
		# total_err -= (target_flux[0] - fluxerr[0]) / (total_comp_flux - total_comp_err)
		# total_err = abs(total_err)
		total_err = 0

		if RANDOM_INSPECTION == True and random.uniform(0, 1) < 0.03:
			plt.imshow(data_minus_bg-loc_tgt_bg, interpolation="nearest", vmin=0.0, vmax=4000.0, cmap="gray", origin="lower")
			plt.colorbar()
			circ = Circle((target_x, target_y), kron_rad, color="b", fill=False)
			ax = plt.gca()
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
		else:
			print("\nNo suitable date header found!!! :(")
			exit()
		MJD.append(mid_photo_time)
		FLUX_RATIO.append(ratio)
		FLUX_RATIO_ERR.append(total_err)
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
		data_file.write("MJD,Flux Ratio\n")
for data_set in COMBINED_DATA:
	MJD.append(data_set[0])
	FLUX_RATIO.append(data_set[1])
	if MAG_MODE == True:
		FLUX_RATIO_ERR.append(0)
	else:
		FLUX_RATIO_ERR.append(data_set[2])
	if EXPORT_CSV == True:
		with open(EXPORT_PATH, "a") as data_file:
			data_file.write("{},{}\n".format(data_set[0], data_set[1]))

the_s = "s" if len(COMPS_SKYCOORD) > 1 else ""

ylabel = f"Flux ratio (target / comp{the_s})" if MAG_MODE == False else "Mangitude"

plt.plot(MJD, FLUX_RATIO, marker="o", color="green", linewidth=0)
# plt.errorbar(MJD, FLUX_RATIO, yerr=FLUX_RATIO_ERR, linewidth=0, elinewidth=2, ecolor="red")
plt.title(f"Time-series flux of {STAR_NAME} in the filter {FILTER}")
plt.xlabel("Time (Days)")
plt.ylabel(ylabel)
plt.show()
