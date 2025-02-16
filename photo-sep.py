# Get a lightcurve from a folder of FITs files

import sep_pjw as sep
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.wcs import WCS
from astropy.time import Time
from os import walk
import math
import numpy as np
import copy
import random
import sys
import time
random.seed(int(time.time()))

# Whether or not to use command line args
USE_CLI_ARGS = True
# other general config
FILE_FOLDER_PATH = "../science-data/"
STAR_NAME = ""
FILTER = "" # filter name as it appears in the FITS headers
MAG_MODE = False # whether or not to use calibrated magnitude instead of instrumental magnitude
RANDOM_INSPECTION = False # randomly inspects a file so you can manually check the apertures
VERBOSE = False # extra logging info in console
EXPORT_FILE = False # the options are: False, "WebObs", or "CSV"
EXPORT_PATH = f"./{STAR_NAME.lower().replace(" ", "-")}-{FILTER.lower()}." + ("csv" if EXPORT_FILE == "CSV" else "txt")
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

# Info of observing site to calculate HJD (lattitude, longitude, altitue)
# If this is not set manually then the code will try to set it automatically from FITS headers
# The "unset" value is -999
SITE_INFO = {
	"lat": -999,
	"lon": -999,
	"alt": -999
	}

if USE_CLI_ARGS:
	TARGET_RA, TARGET_DEC = "", ""
	TARGET_PIX = [-1, -1]
	COMPS_RA, COMPS_DEC = [""], [""]
	COMPS_PIX = [[-1, -1]]
	VERBOSE, MAG_MODE, RANDOM_INSPECTION, EXPORT_FILE = (False,)*4
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
			i -= 1
		elif t_f == "-e": # export data
			EXPORT_FILE = sys.argv[i+1]
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

if MAG_MODE == True and len(COMPS_MAGS) != len(COMPS_RA) and len(COMPS_MAGS) != len(COMPS_PIX):
	raise ValueError("The length of the comparison magnitude array is mismatched.")

if FILTER == "":
	raise ValueError("Filter name has not been set.")

if EXPORT_FILE != False and EXPORT_FILE != "WebObs" and EXPORT_FILE != "CSV":
	raise ValueError("Unknown file export option:", EXPORT_FILE)

if EXPORT_FILE != False and len(COMPS_NAMES) != len(COMPS_RA) and len(COMPS_NAMES) != len(COMPS_PIX):
	raise ValueError("The length of the comparison names array is mismatched.")

if EXPORT_FILE == "WebObs":
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
print("Export photometry data to file:", EXPORT_FILE)
if EXPORT_FILE == "WebObs":
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
			if annulus_mask[y][x] == True and not math.isnan(data[y][x]):
				annulus_pixels.append(data[y][x])
	cutoff = np.percentile(annulus_pixels, 95)
	for i in range(len(annulus_pixels)-1, -1, -1):
		if annulus_pixels[i] > cutoff:
			del annulus_pixels[i]
	return np.mean(annulus_pixels)

# get square slice in 2d array
def getSquareSlice(data, center_x, center_y, side_length):
	data_shape = data.shape
	slice_x = [int(center_x - side_length/2), int(center_x + side_length/2)]
	slice_y = [int(center_y - side_length/2), int(center_y + side_length/2)]

	for i in range(2):
			if slice_x[i] >= data_shape[1]:
				slice_x[i] = data_shape[1] - 1
			elif slice_x[i] < 0:
				slice_x[i] = 0

			if slice_y[i] >= data_shape[0]:
				slice_y[i] = data_shape[0] - 1
			elif slice_y[i] < 0:
				slice_y[i] = 0
	sliced = data[slice_y[0] : slice_y[1], slice_x[0] : slice_x[1]]
	return copy.deepcopy(sliced)

TARGET_SKYCOORD = None
if TARGET_RA != "":
	TARGET_SKYCOORD = SkyCoord(ra=TARGET_RA, dec=TARGET_DEC, unit=(u.hourangle, u.deg))

COMPS_SKYCOORD = []
for i in range(len(COMPS_RA)):
	if COMPS_RA[i] != "":
		COMPS_SKYCOORD.append(SkyCoord(ra=COMPS_RA[i], dec=COMPS_DEC[i], unit=(u.hourangle, u.deg)))

# arrays for exporting later, about target star
MJD = []
BRIGHTNESS = []
BRIGHTNESS_UNCRT = []
fits_filenames = []
for (_, _, filenames) in walk(FILE_FOLDER_PATH):
    fits_filenames.extend(filenames)
    break
fits_filenames = sorted(fits_filenames)

START_TIME = time.time()
TARGET_INST_MAG = [] # target star instrumental mag
COMPS_INST_MAG = [] # first (and second) comparison star instrumental mag
AIRMASS = []

other_star_qtty = max(len(COMPS_SKYCOORD), len(COMPS_PIX)) # comp star plus check star (if present)
comp_qtty = other_star_qtty - 1 if other_star_qtty > 1 and EXPORT_FILE == "WebObs" else other_star_qtty

# add placeholder pixel values to prevent out of bounds access error later
if len(COMPS_PIX) < len(COMPS_SKYCOORD):
	for i in range(0, len(COMPS_SKYCOORD)- len(COMPS_PIX)):
		COMPS_PIX.append([-1, -1])

for file_name in fits_filenames:
	if not file_name.endswith((".fits", ".fits.fz")):
		continue
	file_path = FILE_FOLDER_PATH + file_name
	with fits.open(file_path) as hdul:
		# handle compressed FITS files
		real_hdul = hdul[0] if hdul[0].header["naxis"] != 0 else hdul[1]

		if not "filter" in real_hdul.header:
			continue
		if real_hdul.header["filter"].strip() != FILTER or "-e00." in file_name:
			continue
		print("Opening file path", file_path)

		use_wcs = True
		this_wcs = WCS(real_hdul.header)

		# switch to little endian byte order
		fits_data = real_hdul.data.astype(np.float32)
		if real_hdul.data.dtype.byteorder == ">":
			# numpy < 2.0.0
			# fits_data = real_hdul.data.byteswap().newbyteorder()

			# solution for numpy >= 2.0.0
			fits_data = real_hdul.data.byteswap()
			fits_data = fits_data.view(fits_data.dtype.newbyteorder("<"))

		sep_bkg = sep.Background(fits_data)

		if this_wcs.wcs.ctype[0] != "":
			target_x, target_y = this_wcs.world_to_pixel(TARGET_SKYCOORD)
		else:
			use_wcs = False
			target_x, target_y = TARGET_PIX

		loc_sky_bg = sky_bg_from_annulus(fits_data, target_x, target_y, 20, 5)
		if VERBOSE == True:
			print("Target star\nLocal sky background", loc_sky_bg)
		loc_fits_data = fits_data - loc_sky_bg

		target_x, target_y, _ = sep.winpos(loc_fits_data, [target_x], [target_y], [3.33])
		target_x = target_x[0]
		target_y = target_y[0]
		if VERBOSE == True:
			print("Position: x=", target_x, ", y=", target_y)

		loc_fits_data = getSquareSlice(loc_fits_data, target_x, target_y, 100)

		sources = sep.extract(loc_fits_data, 2.25, err=sep_bkg.globalrms, minarea=15)

		target_i = -1 # index of the target star
		for i in range(len(sources["thresh"])):
			dist = math.sqrt((sources["x"][i] - 50)**2 + (sources["y"][i] - 50)**2)
			if dist < 2:
				target_i = i
				if sources["flag"][i] > 1:
					raise ValueError("SEP raised unexpected flag for target: ", sources["flag"][i])
				break
		if target_i == -1:
			raise ValueError("SEP could not detect target star above threshold!")
		# scale a and b parameters
		target_a = 2.5 * sources["a"][target_i]
		target_b = 2.5 * sources["b"][target_i]
		target_th = sources["theta"][target_i]

		if VERBOSE == True:
			print("Ellipse parameters: a=", target_a, ", b=", target_b, ", theta=", target_th)

		target_flux, target_err, t_flag = sep.sum_ellipse(loc_fits_data, [sources["x"][target_i]], [sources["y"][target_i]],
			[target_a], [target_b], [target_th], err=sep_bkg.globalrms, subpix=10)

		if t_flag[0] != 0:
			raise ValueError("SEP raised a non-zero flag for target star photometry:", t_flag[0])

		target_flux = target_flux[0]
		target_err = target_err[0]

		if VERBOSE == True:
			print("Flux: ", target_flux, "+/-", target_err, "\n")

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
			if VERBOSE == True:
				print("Comparison star", i+1)
				print("Local sky background:", loc_sky_bg)
			loc_fits_data = fits_data - loc_sky_bg
			comp_x, comp_y, cpos_flag = sep.winpos(loc_fits_data, [comp_x], [comp_y], [3])
			comp_x = comp_x[0]
			comp_y = comp_y[0]
			COMPS_PIX[i] = [comp_x, comp_y]
			if VERBOSE == True:
				print("Location: x=", comp_x, ", y=", comp_y)

			loc_fits_data = getSquareSlice(loc_fits_data, comp_x, comp_y, 100)

			sources = sep.extract(loc_fits_data, 2.25, err=sep_bkg.globalrms, minarea=15)

			comp_i = -1 # index of the target star
			for ii in range(len(sources["thresh"])):
				dist = math.sqrt((sources["x"][ii] - 50)**2 + (sources["y"][ii] - 50)**2)
				if dist < 2:
					comp_i = ii
					if sources["flag"][ii] > 1:
						raise ValueError("SEP raised unexpected flag for target: ", sources["flag"][ii])
					break
			if comp_i == -1:
				raise ValueError("SEP could not detect target star above threshold!")

			# scale a and b parameters
			comp_a = 2.5 * sources["a"][comp_i]
			comp_b = 2.5 * sources["b"][comp_i]
			if VERBOSE == True:
				print("Ellipse parameters: a=", comp_x, ", b=", comp_y, ", theta=", sources["theta"][comp_i])

			comp_flux, comp_err, c_flag = sep.sum_ellipse(loc_fits_data, [sources["x"][comp_i]], [sources["y"][comp_i]],
				[comp_a], [comp_b], [sources["theta"][comp_i]], err=sep_bkg.globalrms, subpix=10)
			if c_flag[0] !=0:
				raise ValueError(f"SEP raised a non-zero flag with comp star {i+1} aperture photometry:", c_flag[0])
			comp_flux = comp_flux[0]
			comp_err = comp_err[0]
			if VERBOSE == True:
				print("Flux: ", comp_flux, "+/-", comp_err)

			# calculate instrumental mags for webobs export
			if EXPORT_FILE == "WebObs":
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

			# finally set up a patch for the comp star
			if RANDOM_INSPECTION == True:
				circ2 = Ellipse(xy=(comp_x, comp_y), width=2*comp_a, height=2*comp_b,
					angle=sources["theta"][comp_i]*180./3.1415, color="b", fill=False)
				mat_patches.append(circ2)
		target_uncrt = math.sqrt(target_uncrt)

		if RANDOM_INSPECTION == True and random.uniform(0, 1) < 0.02:
			mean, stdev = np.mean(fits_data), np.std(fits_data)
			fig = plt.figure(figsize=(8, 6))
			ax = fig.add_subplot()
			ax.set_title("Sample image showing target (green) and comparison (blue) stars")
			im1 = ax.imshow(fits_data, interpolation="nearest", vmin=mean-3*stdev, vmax=mean+3*stdev, cmap="gray", origin="lower")
			fig.colorbar(im1)
			circ = Ellipse(xy=(target_x, target_y), width=2*target_a, height=2*target_b, angle=target_th*180./3.1415, color="g", fill=False)
			ax.add_patch(circ)
			for p in mat_patches:
				ax.add_patch(p)
			plt.show()
			RANDOM_INSPECTION = False

		mid_photo_time = 0
		if "DATE-AVG" in real_hdul.header:
			mid_photo_time = Time(real_hdul.header["DATE-AVG"]).jd
		elif "MJD-OBS" in real_hdul.header:
			# mjd is converted to normal jd too
			mid_photo_time = real_hdul.header["MJD-OBS"] + real_hdul.header["EXPTIME"] / 172800 + 2400000.5
		elif "DATE-OBS" in real_hdul.header:
			mid_photo_time = Time(real_hdul.header["DATE-OBS"]).jd + real_hdul.header["EXPTIME"] / 172800
		else:
			raise ValueError("No suitable date header found.")

		if "airmass" in real_hdul.header:
			AIRMASS.append(str(real_hdul.header["airmass"]))
		else:
			AIRMASS.append("na")

		if EXPORT_FILE == "CSV":
			# different variations of the headers
			if "LAT-OBS" in real_hdul.header and SITE_INFO["lat"] == -999:
				SITE_INFO["lat"] = real_hdul.header["LAT-OBS"]
			if "LONG-OBS" in real_hdul.header and SITE_INFO["lon"] == -999:
				SITE_INFO["lon"] = real_hdul.header["LONG-OBS"]
			if "ALT-OBS" in real_hdul.header and SITE_INFO["alt"] == -999:
				SITE_INFO["alt"] = real_hdul.header["ALT-OBS"]

			# another variation is "SITE"
			if "SITELAT" in real_hdul.header and SITE_INFO["lat"] == -999:
				SITE_INFO["lat"] = real_hdul.header["SITELAT"]
			if "SITELONG" in real_hdul.header and SITE_INFO["lon"] == -999:
				SITE_INFO["lon"] = real_hdul.header["SITELONG"]
			if "SITEALT" in real_hdul.header and SITE_INFO["alt"] == -999:
				SITE_INFO["alt"] = real_hdul.header["SITEALT"]
			if "SITEELEV" in real_hdul.header and SITE_INFO["alt"] == -999:
				SITE_INFO["alt"] = real_hdul.header["SITEELEV"]

		MJD.append(mid_photo_time)
		BRIGHTNESS.append(target_brightness)
		BRIGHTNESS_UNCRT.append(target_uncrt)
		TARGET_INST_MAG.append(-2.5 * math.log10(target_flux))
		TARGET_PIX = [target_x, target_y]
		print()

if EXPORT_FILE == "CSV":
	print("Got this observervation site info:", SITE_INFO)

# sort all the data points by MJD so it can be displayed easily
COMBINED_DATA = []
MIN_JD = min(MJD)
for i in range(len(MJD)):
	COMBINED_DATA.append([MJD[i] - MIN_JD,
		BRIGHTNESS[i],
		BRIGHTNESS_UNCRT[i]])
	if EXPORT_FILE == "WebObs":
		idx = len(COMBINED_DATA)-1
		COMBINED_DATA[idx].append(AIRMASS[i])
		COMBINED_DATA[idx].append(TARGET_INST_MAG[i])
		COMBINED_DATA[idx].append(COMPS_INST_MAG[i])

print("Exporting data to file...")

# reset arrays for sorting by MJD
MJD = []
BRIGHTNESS = []
BRIGHTNESS_UNCRT = []
COMBINED_DATA = sorted(COMBINED_DATA, key=(lambda x: x[0]))
# initial webobs file setup
if EXPORT_FILE != False:
	with open(EXPORT_PATH, "a") as data_file:
		file_string = ""
		if EXPORT_FILE == "WebObs":
			file_string = f"#TYPE=EXTENDED\n#OBSCODE={AAVSO_OBSCODE}\n#SOFTWARE=https://github.com/denis-3/SEP-Photometer  (Revision of January 18, 2025)\n"
			file_string += "#DELIM=,\n#DATE=JD\n#OBSTYPE=CCD\n"
			if other_star_qtty == 3:
				file_string += f"#The comparison stars used are {COMPS_NAMES[1]} and {COMPS_NAMES[2]}\n"
			elif other_star_qtty > 3:
				file_string += "#The comparison stars used are " + ", ".join(COMPS_NAMES[1:-1]) + ", and " + COMPS_NAMES[len(COMPS_NAMES)-1] + "\n"

			file_string += "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES"
		elif EXPORT_FILE == "CSV":
			file_string = "JD,HJD," + ("Calibrated" if MAG_MODE == True else "Instrumental") + " Magnitude,Magnitude Uncertainty"
		data_file.write(file_string)


# process data and add it to webobs if necessary
for data_set in COMBINED_DATA:
	MJD.append(data_set[0])
	BRIGHTNESS.append(data_set[1])
	BRIGHTNESS_UNCRT.append(data_set[2])
	if EXPORT_FILE != False:
		with open(EXPORT_PATH, "a") as data_file:
			current_line = ""
			if EXPORT_FILE == "WebObs":
				current_line = f"\n{STAR_NAME},{MIN_JD, data_set[0]},"
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
				current_line += "VMAGINS=" + "%.3f" % data_set[4]
				if other_star_qtty == 2:
					current_line += "|CREFMAG=" + "%.3f" % COMPS_MAGS[1] + "|KREFMAG=" + "%.3f" % COMPS_MAGS[0]
				else:
					current_line += "|KREFMAG=" + "%.3f" % COMPS_MAGS[0]
			elif EXPORT_FILE == "CSV":
				# heliocentric correction
				site_loc = EarthLocation.from_geodetic(lon=SITE_INFO["lon"], lat=SITE_INFO["lat"], height=SITE_INFO["alt"])
				this_time = Time(MIN_JD + data_set[0], format="jd", scale="utc", location=site_loc)
				htime_corr = this_time.light_travel_time(TARGET_SKYCOORD, "heliocentric")

				# JD, HJD, brightness, brightness uncertainty
				current_line += f"\n{this_time},{this_time + htime_corr},{data_set[1]},{data_set[2]}"
			data_file.write(current_line)

print("Displaying chart...")

ylabel = "Instrumental Magnitude" if MAG_MODE == False else "Calibrated Magnitude"

fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot()
ax.set_title(f"Time-series flux of {STAR_NAME} in the filter {FILTER}", fontsize=16)
ax.plot(MJD, BRIGHTNESS, color="green", linewidth=1.5, alpha=0.25)
ax.scatter(MJD, BRIGHTNESS, marker="o", color="green", s=15)
ax.errorbar(MJD, BRIGHTNESS, yerr=BRIGHTNESS_UNCRT, linewidth=0, elinewidth=1.5, ecolor="#333")
ax.set_xlabel(f"Time since {MIN_JD} JD (Days)", fontsize=13)
ax.set_ylabel(ylabel, fontsize=13)
ax.invert_yaxis()

END_TIME = time.time()
TOTAL_TIME = "%.2f" % (END_TIME - START_TIME)
print(f"Complete! That took {TOTAL_TIME} seconds")

plt.show()
