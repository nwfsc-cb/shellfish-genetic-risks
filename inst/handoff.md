# Hand off & thoughts for next steps

I'm handing off the latest version of the model code (``GRs.py``). I've been using Python 3.7.3, and the script dependencies that I had to install (with versions I've been using) are:
- simuPOP v1.1.10.9
- argparse v1.1
- pandas v.0.25.2
- numpy v1.19.2

You can  access the help file with ``python GRs.py -h``, which will return:

```
usage: GRs_chill.py [-h] -b BATCH -r REPS -x COREID [-f FREQ]

An individual-based model of shellfish production, escape, and genetic impacts
to wild populations

optional arguments:
  -h, --help            show this help message and exit
  -b BATCH, --batch BATCH
                        Name of batch to be applied to ouptut files, e.g.,
                        'low_mig'
  -r REPS, --reps REPS  Number of replicates for simulation, e.g., 10
  -x COREID, --coreid COREID
                        Identifier for combining results across cores, e.g., 1
  -f FREQ, --freq FREQ  Path to equilibrium allele frequencies file; if not
                        provided, will initialize biallelic neutral loci to
                        0.5 and 0.5.
```

I'm sure y'all will decide where it needs to go, but I wanted to add my thoughts on next steps for making it more accessible / reproducible:

- Input parameters  
 - Consider inputting parameters values through a parameter input file instead of changing value in the script. Currently most of the parameter values are specified after the functions, other than a few that are taken at the command line ([1] batch name, [2] name for instance and [3] number of replicates in this instance)
  - Clearly link parameter input file to simulation output. Because the parameter values have been specified within the script itself, it was up to my bookkeeping to stay organized. Probably better to hard code that.
  - Consider whether parameters should be a proportion or a value, etc., and making them consistent. For example, there’s stage-based (larvae, juvenile, adults) escape rates and then there’s the number of offspring produced due to gamete escape (different process in the code so I thought about it differently at the time) and that former parameter could be a proportion too.    
- Monthly time steps
 - Currently, there are monthly time steps, yet response variables are measured and reported annually. If you leave as is, you may find mismatch if you manually check things in other months (other than months with reproduction events, which is when simuPOP internally measures the response variables), etc. This also means our response variables don’t always catch seasonal effects. If you wanted you could add code to measure and report the response variables every month, although it could be unnecessary and probably will add a little run time.
 - There are input parameters for the probabilities of different events (e.g., reproduction, escape, seed production) by month, so you can adjust seasonality. Then there are input parameters for the scale of such events, like recruitment size for a reproduction event. If it’s better to have an input parameter that’s the total recruitment per year, and then divide it up seasonally, you could add code for that.    
- Temporal FST. There’s some code that saves the population object at the end of the pre-farm period, and then estimates temporal FST between that population object and the main population object at the end of the simulation. I didn’t end up using this, and you could probably remove this code if no one wants that anymore.

Happy to meet to explain the code!
