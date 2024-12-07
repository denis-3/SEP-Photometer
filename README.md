# SEP Photometer
A high-level photometer program that uses [SEP](https://sep.readthedocs.io/en/latest/). In a word, it iterates over all the FITS files in a given folder to find the target star magnitude or flux ratio over time.

## Installation
Install Python first, then download the `photo-sep.py` file (e.g. with `git clone` or just by copy-pasting). It is recommended to place the file in its own folder.

Next, install the dependencies with this command: `pip install sep astropy matplotlib numpy`

## Usage
First, obtain some FITS images in a single folder, which should be placed preferably in the same directory as the `photo-sep.py` file. It is *very important* that there are *no other files* in the FITS images folder other than the FITS images.

The program can be used with the following command line options. An asterisk **\*** indicates the option is required.

* `--path`**\***: The directory path in which there are FITS images. *Include a trailing slash after the path.*
* `--filter`**\***: The filter of the observation. Only FITS images in the given directory with the value of their "FILTER" header equal to the filter will be used.
* `--tname`**\***: Name of the target star (only used for the title on the graph).
* `--tra`**\***: The Right Ascension of the target star, in `HH MM SS.SS`. The hours, minutes, and seconds components must be separated by a space. E.g. `12 34 56.78`
* `--tdec`**\***: Declination of the target star, in `DD MM SS.SS`. The degree, minutes, and seconds components must be separated by a space. E.g. `-65 43 21.09`
* `--cra`**\***: Right Ascension of the comparison star. It has the same format as `--tra`.
* `--cdec`**\***: Declination of the comparison star. It has the same format as `--tdec`.
* `--cmag`: The magnitude of the comparison star. If this option is omitted, the code will output a flux ratio.
* `-i`: Include this parameter to "inspect" one of the images at random. A graph will show up, displaying the star field and circles around the target and comparison star(s), which can be manually examined for fit. The blue circle represents the aperture around the target star, and the red circle likewise for the comparison star(s).
* `-v`: "Verbose;" include this parameter to output extra information in the console.
* `-e`: Include this option to "export" the data to a CSV file.

Enclose any option values with spaces around quotes `""`, like so: `--tname "Beta Cephei"`.

Here is an example command: `python photo-sep.py --path "./path/to/data/" --filter V --tname "WASP-164b" --tra "22 59 29.67" --tdec "-60 28 48.72" --cra "22 59 48.334" --cdec "-60 30 59.52" -v -e`

After the program finishes analyzing the images, it will display a graph showing the data.

Note: support for multiple comparison stars is available but it must be set up from within the code, not CLI.
