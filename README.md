# A measurement of the Ly-β forest power spectrum and its cross with the Ly-α forest in X-Shooter 100


The Ly-α forest has been used to constrain the Universe's initial conditions and the timing of reionization processes. However, many studies have found that large degeneracies in parameters exist when using the standard statistic, the Ly-α forest power spectrum, such as between the gas temperature and particle mass in warm/fuzzy dark matter models. To break these degeneracies and improve the cosmological constraining power of the forest, we measure the power spectrum of the Ly-β forest and its cross correlation with the coeveal Ly-α forest using the VLT/XSHOOTER XQ-100 Legacy Survey data set, motivated by this transition’s lower absorption cross-section that makes it sensitive to higher densities relative to Ly-α.

# How to run the code (instructions are in progress)

The most important files are `main.py`, `QuasarSpectrum.py` , `inis.py`, `options.py`, and `boot_indo.py`. These are all in the `lyb_pk/pipeline` directory and use data from the `lyb_pk/data` directory. These scripts creates output files with path: `lyb_pk/output/[filename]`. The paths for all outputs are written in the `inis.py` file.

### Easiest test using 100 mocks
1) Make sure you have the directory paths set up the same as mine. `lyb_pk` is home directory and it should contain `data`, `pipeline`, `output`, and probably for your benefit `plot/figures`. Ensure that the mocks are available in the correct filepath. `lyb_pk/data/mocks/XQ-100_lyb_nocorr[filenames]` and that each of the files are zipped. (So they look like this `mock-0000_xq100_lyb_nocorr.txt.gz`)
2) Go to `inis.py` and set `cat_name = "mocks/XQ-100_catalogue_n100"` and `tag=lyb_nocorr`
3) Within the `lyb_pk/pipeline` directory, type into the command line `chmod +x main.py` and then `./main.py` to run the code. 
4) To get error-bars and covariance matrices, run `./boot_indo.py`
5) The final datatable has path `lyb_pk/output/pk_errboot_{0}_{1}.txt".format(mock_or_obs,tag)`
6) If you use `pd.read_csv([file])` on almost any file these scripts create, you will get an easy to read-in table with labeled columns. If you just want the values, you may type, `pd.read_csv([file]).values`
