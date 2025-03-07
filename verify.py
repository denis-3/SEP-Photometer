import matplotlib.pyplot as plt

# settings
KNOWN_FILE_PATH = "./" # file path of past light data
OBS_FILE_PATH = "./" # file path of the transit in question
PERIOD = 1 # period, in days
JD_COLUMN = [1, 1] # Julian date column of data in the first and second file (starting from 1)
BRIGHTNESS_COLUMN = [2, 2] # Brightness column of data in the first and second file (starting from 1)
KNOWN_EPOCH = 1234567.789 # JD of a past known epoch
# the following value is added to the observed brightness.
# useful in cases where known/obs data is in different filters
# and the discrepancy can be corrected by a simple constant
BRIGHTNESS_ADJUST = 0

# K for Known
K_PHASE = []
K_BRIGHTNESS = []
K_TIME = []
# O for Observed
O_PHASE = []
O_BRIGHTNESS = []

files = [open(KNOWN_FILE_PATH, "r"), open(OBS_FILE_PATH, "r")]
for i in range(len(files)):
	for line in files[i]:
		splitted = line.split(",")
		if JD_COLUMN[i] >= len(splitted) or BRIGHTNESS_COLUMN[i] > len(splitted):
			continue
		this_time = None
		this_brightness = None
		try:
			this_time = float(splitted[JD_COLUMN[i] - 1])
			this_brightness = float(splitted[BRIGHTNESS_COLUMN[i] - 1])
		except ValueError:
			continue
		this_phase = ((this_time - KNOWN_EPOCH) / PERIOD) % 1
		if i == 0:
			K_PHASE.append(this_phase)
			K_BRIGHTNESS.append(this_brightness)
			K_TIME.append(this_time)
		elif i == 1:
			O_PHASE.append(this_phase)
			O_BRIGHTNESS.append(this_brightness + BRIGHTNESS_ADJUST)
		else:
			raise ValueError("There are too many files in the array!")


print("Successfully parsed", len(K_PHASE) + len(O_PHASE), "entries")

fig = plt.figure(figsize=(9, 5))
ax = fig.add_subplot()
ax.set_title("Period fold", fontsize=16)
ax.scatter(K_PHASE, K_BRIGHTNESS, marker="o", c=K_TIME, cmap="copper", s=15)
ax.scatter(O_PHASE, O_BRIGHTNESS, marker="o", color="green", s=15)
ax.set_xlabel("Phase", fontsize=13)
ax.set_ylabel("Magnitude", fontsize=13)
ax.invert_yaxis()

plt.show()
