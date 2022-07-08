# MS-Protein-Adduct-Identification
Identification of protein adducts in deconvoluted mass spectrometry data by matching experimental peak distributions with potential isotope patterns.

![GitHub](https://img.shields.io/github/license/dlon450/MS-Binding-Sites-Identification) ![GitHub](https://img.shields.io/badge/Python-v3.9.7-blue)

## Getting Started
To access the web app, run the src/app.py file. 
```
python src/app.py
```
Otherwise, the search can be run directly using src/binding_site_search.py (changing file directories accordingly). 
```
python src/binding_site_search.py
```
Note that **deconvoluted** spectrum will be required as input (for the bound MS csv file).

## Built With
* HTML/JS/CSS
* [Flask](https://flask.palletsprojects.com/en/2.0.x/)
* [PyOpenMS](https://pyopenms.readthedocs.io/en/latest/)
* [OR-Tools](https://developers.google.com/optimization/cp)
* [SimilarityMeasures](https://github.com/cjekel/similarity_measures)

## Contact Information
If you encounter any problems in using the code, please open an issue in this repository.
