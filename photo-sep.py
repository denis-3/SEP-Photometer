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
FILE_FOLDER_PATH = "../science-data/"
STAR_NAME = ""
FILTER = ""
MAG_MODE = False # whether or not to convert flux rations to standardized magnitude
RANDOM_INSPECTION = False # randomly inspects a file so you can manually check the apertures
VERBOSE = False # extra logging info in console
EXPORT_WEBOBS = False # whether or not to export photometry to csv
EXPORT_PATH = f"./{STAR_NAME.lower().replace(" ", "-")}-{FILTER.lower()}.txt"
AAVSO_OBSCODE = "" # observer code, only used in making the WebObs file
AAVSO_FILTER = "" # name of the filter according to the AAVSO abbreviation
AAVSO_CHART = "na" # not applicable by default

# ~~~ Choose a target ~~~ #

TARGET_RA = ""
TARGET_DEC = ""

# if no WCS, insert the initial x, y pixel values of the target star here
TARGET_PIX = [-1, -1]

# ~~~ Choose a comp(s) ~~~ #

# note: in ensemble photometry, the first star in the array is the check star
# note: if only one star is input, there is no check star
COMPS_RA = [""] # string array
COMPS_DEC = [""] # string array
COMPS_MAGS = [0.] # float array
COMPS_NAMES = [""] # string array

# if no WCS, insert the initial x, y pixel values of the comp star(s) here
COMPS_PIX = [[-1, -1]] # array of [int, int]


if USE_CLI_ARGS:
	TARGET_RA, TARGET_DEC = "", ""
	TARGET_PIX = [-1, -1]
	COMPS_RA, COMPS_DEC = [""], [""]
	COMPS_PIX = [[-1, -1]]
	VERBOSE, MAG_MODE, RANDOM_INSPECTION, EXPORT_WEBOBS = (False,)*4
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
			TARGET_PIX = map(int, sys.argv[i+1].split(","))
		elif t_f == "--cra": # comp star RA
			COMPS_RA = sys.argv[i+1].split(",")
		elif t_f == "--cdec": # comp star Dec
			COMPS_DEC = sys.argv[i+1].split(",")
		elif t_f == "--cpix": # comparison star pixel values, separated by space
			pixel_split = sys.argv[i+1].split(",")
			COMPS_PIX = [map(int, item.split(",")) for item in pixel_split]
		elif t_f == "--cname": # comparison star names, separated by comma
			COMPS_NAMES = sys.argv[i+1].split(",")
		elif t_f == "-i": # whether or not to do inspection
			RANDOM_INSPECTION = True
			i -= 1
		elif t_f == "--cmag": # comp star mag
			MAG_MODE = True
			COMPS_MAGS = map(float, sys.argv[i+1].split(","))
		elif t_f == "-v": # verbose
			VERBOSE = True
		elif t_f == "-e": # export data to CSV
			EXPORT_WEBOBS = True
		elif t_f == "--aavso-code": # AAVSO observer code
			AAVOS_OBSCODE = sys.argv[i+1]
		elif t_f == "--aavso-filter-name": # name of the filter, as to report in the AAVSO webobs
			AAVSO_FILTER = sys.argv[i+1]
		elif t_f == "--aavso-chart": # AAVSO chart ID
			AAVSO_CHART = sys.argv[i+1]
		else:
			raise ValueError(f"Unrecrognized option: {t_f}\nPerhaps check for typos in the command and make sure that spaces are enclosed in quotes.")
		i += 2

# various sanity checks on inputs
if TARGET_RA == "" and TARGET_DEC == "" and TARGET_PIX == [-1, -1]:
	raise ValueError("Target star position has not been set.")

if STAR_NAME == "":
	raise ValueError("Target name has not been set.")

if len(COMPS_RA) != len(COMPS_DEC):
	raise ValueError("The length of the comparison RA and Dec arrays are different.")

if COMPS_RA[0] == "" and COMPS_PIX[0] == [-1, -1]:
	raise ValueError("Both the comparison celestial and pixel coordinates were not input.")

if MAG_MODE == True and len(COMPS_MAGS) != len(COMPS_RA):
	raise ValueError("The length of the comparison magnitude array is mismatched.")

if FILTER == "":
	raise ValueError("Filter name has not been set.")

if EXPORT_WEBOBS == True and len(COMPS_NAMES) != len(COMPS_RA):
	raise ValueError("The length of the comparison names array is mismatched.")

if EXPORT_WEBOBS == True:
	if AAVSO_OBSCODE == "":
		raise ValueError("No AAVSO Observer Code has been input, while the 'export to WebObs' option is on")

	if AAVSO_FILTER == "":
		raise ValueError("No AAVSO filter abbreviation has been provided, while the 'export to WebObs' option is on")

	if AAVSO_CHART == "na":
		print("\033[1;33;40m~~~\nThe AAVSO chart ID is set to 'na'. It is recommended to input an AAVSO chart ID corresponding to the starfield.\n~~~\033[0m")


print("Initializing SEP photometer with these settings:")
print("File folder path:", FILE_FOLDER_PATH)
print("Filter:", FILTER)
print("Target star RA:", TARGET_RA)
print("Target star Dec:", TARGET_DEC)
print("Comp star(s) RA:", COMPS_RA)
print("Comp star(s) Dec:", COMPS_DEC)
print("Comp star pixel coordinates:", COMPS_PIX)
print("Comp star name(s):", COMPS_NAMES)
print("Inspection:", RANDOM_INSPECTION)
print("Magnitude mode:", MAG_MODE)
if MAG_MODE == True:
	print("Comp star(s) Magnitudes", COMPS_MAGS)
print("Export to WebObs Extended File Format:", EXPORT_WEBOBS)
if EXPORT_WEBOBS == True:
	print("WebObs export path:", EXPORT_PATH)
	print("AAVSO Observer Code:", AAVSO_OBSCODE)
	print("AAVSO filter identifier:", AAVSO_FILTER)
	print("AAVSO Star Chart ID:", AAVSO_CHART)
	print("Note: the first \"comp\" star is actually a check star")
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
			if annulus_mask[y][x] == True:
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

# arrays for exporting later, about target star
MJD = []
BRIGHTNESS = []
BRIGHTNESS_UNCRT = []

fits_filenames = []
for (_, _, filenames) in walk(FILE_FOLDER_PATH):
    fits_filenames.extend(filenames)
    break

START_TIME = time.time()
TARGET_INST_MAG = [] # target star instrumental mag
COMPS_INST_MAG = [] # first (and second) comparison star instrumental mag
AIRMASS = []

other_star_qtty = max(len(COMPS_SKYCOORD), len(COMPS_PIX)) # comp star plus check star (if present)
comp_qtty = other_star_qtty - 1 if other_star_qtty > 1 and EXPORT_WEBOBS == True else other_star_qtty

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
		# numpy < 2.0.0
		# fits_data = real_hdul.data.byteswap().newbyteorder()

		# solution for numpy >= 2.0.0
		fits_data = real_hdul.data.byteswap()
		fits_data = fits_data.view(real_hdul.data.dtype.newbyteorder("<"))


		if file_name.endswith(".fits.fz"):
			fits_data = real_hdul.data
		sep_bkg = sep.Background(fits_data)
		sep_bkg_err = sep_bkg.globalrms

		if this_wcs.wcs.ctype[0] != "":
			target_x, target_y = this_wcs.world_to_pixel(TARGET_SKYCOORD)
		else:
			use_wcs = False
			target_x, target_y = TARGET_PIX

		loc_sky_bg = sky_bg_from_annulus(fits_data, target_x, target_y, 12, 5)
		loc_fits_data = fits_data - loc_sky_bg
		target_x, target_y, _ = sep.winpos(loc_fits_data, [target_x], [target_y], [3.33])
		target_x = target_x[0]
		target_y = target_y[0]

		# automatic radius - captures most of the flux of the star
		auto_rad = sep.flux_radius(loc_fits_data, [target_x], [target_y], [12], 0.9)[0][0]
		target_flux, target_err, t_flag = sep.sum_circle(loc_fits_data, [target_x], [target_y], [auto_rad],
			err=sep_bkg_err, subpix=10)

		if t_flag[0] != 0:
			raise ValueError("\nSEP raised a non-zero flag for target star photometry:", t_flag[0])

		target_flux = target_flux[0]
		target_err = target_err[0]

		if VERBOSE == True:
			print("Target X Y coords and radius", target_x, target_y, auto_rad)

		target_brightness = 0 # target brightness
		target_uncrt = 0 # target uncertainty
		mat_patches = [] # matplotlib patches
		for i in range(max(len(COMPS_SKYCOORD), len(COMPS_PIX))):
			if use_wcs:
				comp_coord = COMPS_SKYCOORD[i]
				comp_x, comp_y = this_wcs.world_to_pixel(comp_coord)
			else:
				comp_x, comp_y = COMPS_PIX[i]

			loc_sky_bg = sky_bg_from_annulus(fits_data, comp_x, comp_y, 12, 5)
			loc_fits_data = fits_data - loc_sky_bg
			comp_x, comp_y, _ = sep.winpos(loc_fits_data, [comp_x], [comp_y], [3.33])
			comp_x = comp_x[0]
			comp_y = comp_y[0]

			comp_arad = sep.flux_radius(loc_fits_data, [comp_x], [comp_y], [12], 0.9)[0][0]
			comp_flux, comp_err, c_flag = sep.sum_circle(loc_fits_data, [comp_x], [comp_y], [auto_rad],
				err=sep_bkg_err, subpix=10)

			if c_flag[0] !=0:
				raise ValueError(f"\nSEP raised a non-zero flag with comp star {i+1} aperture photometry:", c_flag[0])

			comp_flux = comp_flux[0]
			comp_err = comp_err[0]

			if VERBOSE == True:
				print("Comp X Y coords and radius", comp_x, comp_y, comp_arad)

			# calculate instrumental mags for webobs export
			if EXPORT_WEBOBS == True:
				if i == 0:
					COMPS_INST_MAG.append([-2.5 * math.log10(comp_flux)])
					continue # check star is not inlcuded in comp ensemble
				elif i == 1 and other_star_qtty == 2:
					COMPS_INST_MAG[len(COMPS_INST_MAG)-1].append(-2.5 * math.log10(comp_flux))

			tc_ratio = target_flux / comp_flux # target-comp ratio
			tc_min_ratio = (target_flux - target_err) / (comp_flux + comp_err)
			tc_max_ratio = (target_flux + target_err) / (comp_flux - comp_err)

			# calculate instrumental magnitudes first
			this_target_mag = - 2.5 * math.log10(tc_ratio)
			min_target_mag = - 2.5 * math.log10(tc_max_ratio)
			max_target_mag = - 2.5 * math.log10(tc_min_ratio)
			target_mag_uncrt = max(max_target_mag - this_target_mag, this_target_mag - min_target_mag)

			if MAG_MODE == True:
				this_target_mag += COMPS_MAGS[i]

			target_brightness += this_target_mag / comp_qtty
			target_uncrt += target_mag_uncrt ** 2 # uncertainty is summed in quadrature

			COMPS_PIX[i] = [comp_x, comp_y]

			# finally set up a patch for the comp star
			if RANDOM_INSPECTION == True:
				circ2 = Circle((comp_x, comp_y), comp_arad, color="r", fill=False)
				mat_patches.append(circ2)
		target_uncrt = math.sqrt(target_uncrt)

		if RANDOM_INSPECTION == True and random.uniform(0, 1) < 0.02:
			fig = plt.figure(figsize=(8, 6))
			ax = fig.add_subplot()
			ax.set_title(f"Sample image of target star (blue) and comps (red)")
			im1 = ax.imshow(loc_fits_data, interpolation="nearest", vmin=0.0, vmax=4000.0, cmap="gray", origin="lower")
			fig.colorbar(im1)
			circ = Circle((target_x, target_y), auto_rad, color="b", fill=False)
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
			raise ValueError("No suitable date header found.")

		if "airmass" in real_hdul.header:
			AIRMASS.append(str(real_hdul.header["airmass"]))
		else:
			AIRMASS.append("na")
		MJD.append(mid_photo_time)
		BRIGHTNESS.append(target_brightness)
		BRIGHTNESS_UNCRT.append(target_uncrt)
		TARGET_INST_MAG.append(-2.5 * math.log10(target_flux))
		TARGET_PIX = [target_x, target_y]
		print()

# sort all the data points by MJD so it can be displayed easily
COMBINED_DATA = []
MIN_MJD = min(MJD)
for i in range(len(MJD)):
	COMBINED_DATA.append([MJD[i] - MIN_MJD,
		BRIGHTNESS[i],
		BRIGHTNESS_UNCRT[i]])
	if EXPORT_WEBOBS == True:
		idx = len(COMBINED_DATA)-1
		COMBINED_DATA[idx].append(AIRMASS[i])
		COMBINED_DATA[idx].append(TARGET_INST_MAG[i])
		COMBINED_DATA[idx].append(COMPS_INST_MAG[i])


# reset arrays for sorting by MJD
MJD = []
BRIGHTNESS = []
BRIGHTNESS_UNCRT = []
COMBINED_DATA = sorted(COMBINED_DATA, key=(lambda x: x[0]))
# initial webobs file setup
if EXPORT_WEBOBS == True:
	with open(EXPORT_PATH, "a") as data_file:
		file_string = f"#TYPE=EXTENDED\n#OBSCODE={AAVSO_OBSCODE}\n#SOFTWARE=https://github.com/denis-3/SEP-Photometer  (Revision of January 17, 2025)\n"
		file_string += "#DELIM=,\n#DATE=JD\n#OBSTYPE=CCD\n"
		if other_star_qtty == 3:
			file_string += f"#The comparison stars used are {COMPS_NAMES[1]} and {COMPS_NAMES[2]}\n"
		elif other_star_qtty > 3:
			file_string += f"#The comparison stars used are {", ".join(COMPS_NAMES[1:-1])}, and {COMPS_NAMES[2]}\n"

		file_string += "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES"
		data_file.write(file_string)

# process data and add it to webobs if necessary
for data_set in COMBINED_DATA:
	MJD.append(data_set[0])
	BRIGHTNESS.append(data_set[1])
	BRIGHTNESS_UNCRT.append(data_set[2])
	if EXPORT_WEBOBS == True:
		with open(EXPORT_PATH, "a") as data_file:
			current_line = f"\n{STAR_NAME},{data_set[0] + MIN_MJD},"
			mag_to_write = data_set[1]
			mag_err_to_write = data_set[2]
			# if standardized magnitudes are not available, we use differential magnitude
			if MAG_MODE == False:
				mag_to_write = -2.5 * math.log10(data_set[1])
				mag_err_to_write = 2.5 / math.log(10) * data_set[2]
			current_line += "%.3f" % mag_to_write
			current_line += ","
			current_line += "%.6f" % mag_err_to_write
			current_line += f",{AAVSO_FILTER},NO,"
			if MAG_MODE == True:
				current_line += "STD"
			else:
				current_line += "DIFF"

			# comps and check star listing
			if len(COMPS_RA) == 1:
				current_line += f",{COMPS_NAMES[0]},{data_set[5][0]},na,na"
			elif len(COMPS_RA) == 2:
				current_line += f",{COMPS_NAMES[1]},{data_set[5][1]},{COMPS_NAMES[0]},{data_set[5][0]}"
			elif len(COMPS_RA) > 2:
				current_line += f",ENSEMBLE,na,{COMPS_NAMES[0]},{data_set[5][0]}"
			current_line += f",{data_set[3]},na,{AAVSO_CHART}," # airmass, group, AAVSO chart ID
			# add notes now
			current_line += f"VMAGINS={"%.3f" % data_set[4]}"
			if other_star_qtty == 2:
				current_line += f"|CREFMAG={"%.3f" % COMPS_MAGS[1]}|KREFMAG={"%.3f" % COMPS_MAGS[0]}"
			else:
				current_line += f"|KREFMAG={"%.3f" % COMPS_MAGS[0]}"

			data_file.write(current_line)


ylabel = "Instrumental Magnitude" if MAG_MODE == False else "Calibrated Magnitude"

fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot()
ax.set_title(f"Time-series flux of {STAR_NAME} in the filter {FILTER}")
ax.plot(MJD, BRIGHTNESS, color="green", linewidth=1.5, alpha=0.25)
ax.scatter(MJD, BRIGHTNESS, marker="o", color="green", s=15)
ax.errorbar(MJD, BRIGHTNESS, yerr=BRIGHTNESS_UNCRT, linewidth=0, elinewidth=1.5, ecolor="#333")
ax.set_xlabel("Time (Days)")
ax.set_ylabel(ylabel)
ax.invert_yaxis()

END_TIME = time.time()
TOTAL_TIME = "%.2f" % (END_TIME - START_TIME)
print(f"Complete! That took {TOTAL_TIME} seconds")

plt.show()
