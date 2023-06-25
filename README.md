# ribozero_Plus_module

* needs a working conda installation and conda in path

* Snakemake needs to be installed in the working enviroment

* sample_working.tsv provides a example of the needed sampleseheet

* use the module like this:

```
$ cd ribozero_Plus_module
```

the current working dirtectory needs to be the rna_module_1 folder

```
snakemake --use-conda --cores 8 --conda-prefix "folder where conda should store the envs data (is optional)"

```


* the pipeline creates a folder structure and writes all the output in the parent directory of the current working directory
