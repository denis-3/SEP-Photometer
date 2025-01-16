
# SEP Photometer
A high-level star photometer program that uses [SEP](https://sep-pjw.readthedocs.io/en/latest/). In a word, it iterates over all the FITS files in a given folder to find the target star magnitude or flux ratio over time. It can also export data to an image file and the [WebObs format](https://www.aavso.org/aavso-extended-file-format).

## Installation
Install Python first, then download the `photo-sep.py` file (e.g. with `git clone` or just by copy-pasting). It is recommended to place the file in its own folder.

Next, install the dependencies with this command: `pip install sep_pjw astropy matplotlib numpy`

## Usage
First, obtain some FITS images in a single folder, which should be placed preferably in the same directory as the `photo-sep.py` file. It is *very important* that there are *no other files* in the FITS images folder other than the FITS images.

The program can be used with the following command line options. An asterisk **\*** indicates the option is required.

* `--path`**\***: The directory path in which there are FITS images. *Include a trailing slash after the path*, e.g. `/home/myuser/Documents/data/`.
* `--filter`**\***: The filter of the observation. Only FITS images in the given directory with the value of their "FILTER" header equal to the filter will be used.
* `--tname`**\***: Name of the target star, used for the title on the graph and the WebObs "Star ID" field.
* `--tra`**\***: The Right Ascension of the target star, in `HH MM SS.SS`. The hours, minutes, and seconds components must be separated by a space. E.g. `12 34 56.78`.
* `--tdec`**\***: Declination of the target star, in `DD MM SS.SS`. The degree, minutes, and seconds components must be separated by a space. E.g. `-65 43 21.09`.
* `--tpix`**\***: The x and y pixel values of the target star, if no WCS is available. It should be input as a string with the format `"x,y"`, e.g. `"100,200"`. Only required if the `--tra` and `--tdec` parameters are not set.
* `--cra`**\***: Right Ascension of the comparison star. It has the same format as `--tra`. It may also be input as a comma-separated string, e.g. `"12 34 56.78,24 68 02.46"`.
* `--cdec`**\***: Declination of the comparison star. It has the same format as `--tdec`. It may also be input as a comma-separated string, e.g. `"-65 43 21.09,08 64 20.86"`.
* `--cpix`**\***: The x and y pixel values of the comparison star, if no WCS is available. Only required if the `--cra` and `--cdec` parameters are not set. It should be input as a string with the format `"x1,y1 x2,y2 x3,y3 ..."` (*note the spaces*), e.g. `"200,400 600, 800 1000,1200"`.
* `--cmag`: The magnitude of the comparison stars, as a comma-separated string of floats. If this option is omitted, the code will output a flux ratio.
* `-i`: Include this parameter to "inspect" one of the images at random. A graph will show up, displaying the star field and circles around the target and comparison star(s), which can be manually examined for fit. The blue circle represents the aperture around the target star, and the red circle likewise for the comparison star(s).
* `-v`: "Verbose;" include this parameter to output extra information in the console.
* `-e`: Include this option to "export" the data to a WebObs file.
* `--aavso-code`: The AAVSO Observer Code of the scientist. Check this page for more information: https://dev-mintaka.aavso.org/new-observer-FAQs . *It is required to set this parameter when using `-e`*.
* `--aavso-filter`: The value of the "Filter" field in the [WebObs](https://www.aavso.org/aavso-extended-file-format) file. *It is required to set this parameter when using `-e`*.
* `--aavso-chart`: The AAVSO chart ID. This parameter is not required to set, but encouraged.

Enclose any option values with spaces around quotes `""`, like so: `--tname "Beta Cephei"`.

Here is an example command: `python photo-sep.py --path "./path/to/data/" --filter V --tname "WASP-164b" --tra "22 59 29.67" --tdec "-60 28 48.72" --cra "22 59 48.334" --cdec "-60 30 59.52" -v -e`

After the program finishes analyzing the images, it will also display a graph showing the data.
