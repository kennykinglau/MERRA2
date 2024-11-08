# Installing am Atmospheric Model
Follow [these instructions](https://lweb.cfa.harvard.edu/~spaine/am/download/src/INSTALLING).

# Running `merra2Player`
1. Create a [NASA Earthdata](https://urs.earthdata.nasa.gov/home) account.
> [!WARNING]  
> Your password must be stored in plaintext on the computer running
> `merra2Player`. Therefore, it is strongly recommended that you choose a new
> password not used on any other website.
2. Approve the GES DISC application for authorization.
   <ol>
     <li>Navigate to Applications > Authorized Apps</li>
     <li>Click on "Approve More Applications"</li>
     <li>Search for "NASA GESDISC DATA ARCHIVE" (without the quotes)</li>
     <li>Click on "Authorize" next to "NASA GESDISC DATA ARCHIVE" and follow the steps.</li>
   </ol>
3. Install Wget (if not already installed)
> [!TIP]
> Wget is available with most package managers. If you're using macOS, it is
> recommended to first install [Homebrew](https://brew.sh/) and then install
> [Wget from Homebrew](https://formulae.brew.sh/formula/wget). You can find
> other installation methods
> [here](http://wget.addictivecode.org/FrequentlyAskedQuestions.html#download).
4. Follow [these steps](https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Generate%20Earthdata%20Prerequisite%20Files)
   to generate Earthdata prerequisite files.
5. [Create](https://docs.python.org/3/library/venv.html#creating-virtual-environments) and [activate](https://docs.python.org/3/library/venv.html#how-venvs-work) a virtual environment with venv.
6. Install requirements with `pip install -r merra2_requirements.txt`
7. Set up directories
> [!NOTE]
> This step should be done from the `python` directory and may not be necessary
> from within the Harvard CANNON environment.
```sh
mkdir -p kovac_lab/keck/wvr_products/merra2_analysis
ln -s kovac_lab/keck/wvr_products/merra2_analysis merra2_products
mkdir merra2_products/merra2_raw_data
mkdir merra2_products/tipper_raw_data
ln -s <path to am executable> am
```

You're now ready to run merra2Player, which can be done in a few ways.

- `predict_Tsky.py`

`predict_Tsky.py` is a command-line tool that drives `merra2Player`. It must be
run from within the `python` directory.
```sh
# Get help on the options available
./predict_Tsky.py --help
# Run an example
./predict_Tsky.py -l SouthPole -d "20230501, 20230503"
```
You can replace SouthPole in the above example with any of the following named
sites (specified in `merra2Player.defineSite`):
- SouthPole
- ChajnantorPlateau
- ChajnantorCerro
- MaunaKea
- Summit
- Qubic

`predict_Tsky.py` can also optionally take a comma-separated list for the `-l`
option specifying any set of coordinates. However, certain atmospheric model
parameters (e.g. volume mixing ratio) are interpolated based on their values at
named sites when coordinates are specified in this way, and the interpolation
can result in impossible values for some parameters (e.g. a volume mixing ratio
less than 0 or greater than 1), which will cause `am` to fail. This failure
is more likely the further away the coordinates are from one of the named sites.

```sh
# Succeeds as of November, 2024. Interpolation works due to proximity to SouthPole.
./predict_Tsky.py -l SouthPoleCoordinates,SouthPoleCoordinates,-89,0,2835 -d "20230501, 20230503"
# Fails as of November, 2024 (volume mixing ratio out of range). Not close to any named sites.
./predict_Tsky.py -l KittPeak,KittPeak,31.9583,-111.5967,2096 -d "20230501, 20230503"
```

- Python interpreter

Start a Python interpreter from within the `python` directory. You should be
able to import the class and call its methods directly, like in the below
example.
```python
import merra2Player as m2p
m = m2p.merra2Player()
m.runMERRA(dateopt={'start': '20230501', 'end': '20230503'})
```
