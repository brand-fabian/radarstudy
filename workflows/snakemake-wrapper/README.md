# snakemake-wrappers

In-house snakemake wrapper repository.

## Usage

To use the wrappers, `clone` this repository and then use the following
statement in your snakemake rules instead of `shell` or `script`:

    wrapper: "file:relative/path/to/repository/and/wrapper"

or

    wrapper: "file:///path/to/wrapper"

## Contributing

Please feel free to contribute to this repository. Your contribution
could either be a new wrapper or unit tests for existing wrappers. Keep 
in mind that we ideally want to synchronize these wrappers with the 
existing wrappers ([here](https://snakemake-wrappers.readthedocs.io/en/stable/))
from other people.