# Meteosat Third Generation

The MTG satellite provides input data for our code. This is a geostationary satellite scanning the globe every 15 minutes.

By default, specMAGIC is using the 640 nm band. The code can be run with other bands (in MTG Normal Resolution, i.e. 1km x 1km), but this is at the user's own risk. Validation of results with these different bands is still ongoing work.

## Resources

A huge amount of useful information on using the MTG data, what it is and how to interpret it, [is provided here](https://user.eumetsat.int/resources/user-guides/mtg-fci-level-1c-data-guide) in the user guide.

There exists [a video recording ](https://www.youtube.com/watch?v=msLlQDhEvZY)of a webinar on the FCI level 1C data products.

Satpy [documentation](https://satpy.readthedocs.io/en/latest/examples/fci_l1c_natural_color.html) for generating a natural colour RBG image from MTG data.

A [more detailed overview](https://resources.eumetrain.org/data/7/739/slides/20250626-HOLL-presentation-data-visualisation-pytroll-satpy.pdf) of using Satpy/pytroll for visualisation/processing of MTG data.

A [Jupyter notebook ](https://gitlab.eumetsat.int/eumetlab/data-services/eumdac_data_store/-/blob/master/1_5_MTG_FCI_data_access.ipynb)on downloading MTG data from the EUMETSAT data store and visualising it.