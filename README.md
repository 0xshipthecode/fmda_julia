fmda_julia
==========

This code is part of a larger project focused on data assimilation in fire weather.

The purpose of this code is to take the output of the Weather Research and Forecasting code [WRF](http://www.wrf-model.org/ WRF)
together with remote automatic weather station observations (RAWS) and produce improved estimates of dead fuel moisture in the
area covered by the WRF grid.

The fmda_julia code uses a Trend Surface Modeling approach together with a Kalman filter running at each grid point to
assimilate ground station observations from the [MesoWest](http://http://mesowest.utah.edu/ "MesoWest") network into the moisture model.

The code simulates three types of dead fuel: 1hr, 10hr and 100hr.

Prerequisites
-------------

  * [julia](http://julialang.org "julia") v0.2.0 (works with `0.2.0-2087.r8e1e64cb`) with the `Calendar` package
  * [fmda_scraper](https://github.com/vejmelkam/fmda_scraper "fmda_scraper")
  * a `wrfout` file from a simulation, using either [WRF](http://www.wrf-model.org/index.php "WRF") code or [WRF-fire](http://openwfm.org "WRF-fire") code (recommended)

**NOTE**: [julia](http://julialang.org "julia") is a very young language that is evolving fast, so using a different version than above may break the code.

Installing Prerequisites
------------------------

**Julia**

    git clone http://github.com/JuliaLang/julia.git
    cd julia
    git checkout ded337b9
    make

The make step will take some time as many dependencies are pulled in and compiled.  After the compilation, the ''julia'' executable file should be available in the current directory.

**fmda_scraper**

Download [fmda_scraper](https://github.com/vejmelkam/fmda_scraper "fmda_scraper").  No installation is needed, the scraper can be used from the git directory.


Obtaining MesoWest RAWS observations
------------------------------------

The user must have a list of station ids which are pertinent to the given region before using the scraper or the assimilation code.

RAWS measurements for the US can be obtained from the [MesoWest](http://mesowest.utah.edu/ MesoWest) website and must be
converted to the format consumed by this code.  This task can be accomplished using the [fmda_scraper](http://github.com/vejmelka/fmda_julia "fmda_julia") code, [documentation](https://github.com/vejmelkam/fmda_scraper/blob/master/README.md "documentation").

The station `info` and `obs` files should be available on the filesystem where the fmda_julia code should run.


Setting up a configuration
--------------------------

The fmda_julia requires a configuration file to run.  An example configuration file is:

    [
        "station_info_dir" => "../station_infos",
        "station_data_dir" => "../station_data",
        "station_info" => "../station_infos/station_list",
        "output_dir" => "test_output",
        "wrf_output" => "../wrf/WRFV3/run/wrfout_d01_2013-06-04_00:00:00",
        "Q" => [1e-4, 5e-5, 1e-5, 0, 0, 0, 5e-5, 5e-5, 0],
        "P0" => [0.01, 0.01, 0.01, 0, 0, 0, 0.01, 0.01, 0],
        "covariates" => [ :constant, :temperature, :pressure, :rain, :lon, :lat, :elevation ],
        "assimilation_time_window" => 3600   # assimilate everything in a 1hr window around the model time
    ]

This file is directly readable by julia and parsing it results in a dictionary that is queried for configuration options
at run-time.

The first two keys `station_info_dir`, `station_data_dir` point to the directories, where the `info` and `obs` files
are stored.  The third key `station_info` is the list of stations that should be used by the code.  The format of the
list is identical to that for [fmda_scraper](http://github.com/vejmelkam/fmda_scraper "fmda_scraper").  The `output_dir`
points to a directory where the execution log and output files will be stored.

The ''output_dir'' must exist.

The ''wrf_output'' file is the WRF output file for the domain on which the moisture code is to be run.  **NOTE** always have a backup of the wrf_output file or use a copy of the output as the wrfout file is modified by the moisture code.



Model settings
--------------

The `Q` parameter is the diagonal of the process noise matrix in the Kalman filters.  The process noise matrix is a
diagonal 9x9 matrix, which represents the state covariance increase incurred by running the moisture model one step.
It is recommended to leave zero values in the vector and modify the rest of the values according to the length of the
time step.  A discussion of process noise is beyond the scope of this documentation.

The `P0` parameter is a the diagonal of a diagonal matrix representing the uncertainty of the background guess.  In the
fmda_julia code, the first estimate of the moisture state is the atmospheric current moisture equilibrium. This background
guess should carry with it a large uncertainty as it is likely to be incorrect.

The `covariates` parameter represents the covariates used in the linear regression part of the trend surface model.  There
is no reason to change this unless very few observation stations are available for the simulated area.  The number of covariates should be less than or equal to the number of obsevations available at each time step.  In the event tht the number of
covariates is equal to the number of observations, the trend surface model will be able to fit the observations exactly.

The `assimilation_time_window` defines the validity of observations in seconds.  Each observation is valid for the duration 
of the `assimilation_time_window` centered on the timestamp of the observation.

Running the julia code
----------------------

**IMPORTANT** Always use a copy of the wrfout file to run the fuel moisture assimilation code.  The code modifies the file itself and thus there exists a danger that the file may be corrupted (if the moisture code crashes).


Execute the following in the `fmda_julia` directory:

    ./run_data_assimilation.jl <path_to_config_file>

This will load all available inputs and start processing the wrfout file and the obtained measurements.

Results of the run
------------------

The julia code will create three new variables in the NetCDF wrfout file, named ''FM1'', ''FM10'' and ''FM100'', which will contain the assimilated fuel moisture.  The first time instant is empty, the fuel moisture values are filled from the second time instant forward.
