# MS-Binding-Sites-Identification
Identification of protein binding sites in deconvoluted mass spectrometry data by matching experimental peak distributions with potential isotope patterns.

![GitHub](https://img.shields.io/github/license/dlon450/MS-Binding-Sites-Identification)

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

## Contact Information
If you encounter any problems in using the code, please open an issue in this repository.
