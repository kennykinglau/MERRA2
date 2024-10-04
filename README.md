# Installing am Atmospheric Model
Follow [these instructions](https://lweb.cfa.harvard.edu/~spaine/am/download/src/INSTALLING).

# Running merra2Player
1. Create a [NASA Earthdata](https://urs.earthdata.nasa.gov/home) account.
> [!WARNING]  
> Your password must be stored in plaintext on the computer running
> merra2Player. Therefore, it is strongly recommended that you choose a new
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
mkdir -p kovac_lab/keck/wvr_products/merra2_analysis/SouthPole
ln -s kovac_lab/keck/wvr_products/merra2_analysis merra2_products
ln -s <path to am executable> am
```
8. Run class from Python interpreter inside the python directory
```python
import merra2Player as m2p
m = m2p.merra2Player()
m.runMERRA(dateopt={'start': '20230501', 'end': '20230503'})
```
