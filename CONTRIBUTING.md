# Contributing 

We welcome all contributions to GeoBO. Contributions can be in the form of new code functions, improvements to the documentation, or by pointing out a bug or potential improvement.  


### Bugs reports and feature requests

For bugs and suggestions, please submit a bug report or feature request use the [Github issue tracker](https://github.com/sebhaan/geobo/issues). Search for existing and closed issues first. If your problem or idea is not yet addressed, please open a new issue. Github
allows you to classify your issues so that we know if it is a bug report, feature request or feedback to the authors.


### Submitting code

If you would like to contribute, please submit a pull request. (See the [Github Hello World](https://guides.github.com/activities/hello-world/) example, if you are new to Github). Please document your functions clearly.

By contributing to the repository you state you own the copyright to those contributions and agree to include your contributions as part of this project under the APGL license.

## Future Work

We welcome all contributions to GeoBO, in particular we would appreciate all pull requests pertaining to the following areas:

- Implementing new forward models in `sensormodel.py`
- Expanding the example studies for other use cases


### Testing

(Currently in development).
Continuous integration and testing is performed with `pytest`, which means running this in the source directory:

```bash
python setup.py test
```

The existing tests should be passing before you start coding (raise an issue if that is not the case) and any new functions should also have tests that we can use to verify the code. 

