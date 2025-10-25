# BlasterSim

BlasterSim simulates pneumatic and spring compressed gas guns.

## Installation

Note: At present, no binary or source distributions are available yet. The following sections were written in anticipation of the future release of the BlasterSim beta version.

## From the binary

The `blastersim` executable on Linux and the `blastersim.exe` executable on Windows can either be placed in the directory you want to run BlasterSim from or placed on your `PATH`.

## From source

The provided source distribution has no dependencies other than a modern Fortran compiler and a Make program like GNU Make on Linux or [jom](https://wiki.qt.io/Jom) or [NMAKE](https://learn.microsoft.com/en-us/cpp/build/reference/nmake-reference?view=msvc-170) on Windows. Unfortunately, BlasterSim will not compile with every Fortran compiler due to the units module being too large. BlasterSim will compile with gfortran and for that reason gfortran is recommended.

After extracting the provided source archive, `cd` into the directory and build BlasterSim. On Linux you can type:

    make blastersim

On Windows you can type (if using jom):

    jom blastersim.exe

Now a BlasterSim executable has been built and can be placed appropriately as discussed in the "From the binary" section.

## From Git

Building from the Git repository requires first building genunits from Ben Trettel's [FLT](https://github.com/btrettel/flt/) repository. The requirement for a modern Fortran compiler remains the same as that from the 

Clone the FLT repository and build genunits:

    git clone https://github.com/btrettel/flt.git
    cd flt
    make genunits

On Windows with jom, you should replace `make genunits` with `jom genunits.exe`.

Place `genunits` or `genunits.exe` on your `PATH`. Now you can build BlasterSim by cloning the BlasterSim repository, `cd`ing into the clone, and following the instructions in the "From source" section.

## Usage

TODO
