# beetlesNMR

Analysis tools for NMR analysis from the univeristy of Liverpool. this package was developed for my own NMR metabolomics analysis and other NMR users might find this package helpful in their future data analysis as it streamlines the data analysis workflow. please report any issues or suggestions which you think might be add to usefulness of the package

### Installation
to instlall package copy and paste the following to your R console
```
devtools::install_github('mikies21/beetlesNMR')
```

### Reading in NMR data

in my previous analysis I have used a tool insude the galaxy.ac.uk universe which allows users, by zipping the bruker experiment files and uploading them on the platform, to extract the complete raw spectra of all the experiments (ppm in the first column and samplesID in the subsequent columns).
To automate this a function was created (**however not completely optimised yet. reads the original raw data and any changes (ie. alignment) made with topspin will not be accounted. need help fixing it plssss**)

```
raw_NMR_df <- NMRMetab_readBruker('path_to experiment_folder')
```
