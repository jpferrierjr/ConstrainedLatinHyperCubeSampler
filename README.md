# Constrained Latin HyperCube Sampler

This version allows for a constrained relationship between two variables. I want to expand this to other variables but it's currently only able to handle two.

Example of how to use the class:

```python
    samples         = 25                            # The amount of LHCS samples requested

    Temperature     = np.linspace( 300, 1100, 50, dtype = np.float32 )
    Pressure        = 0.0011875*( Temperature - 300 )**2

    datanames       = ['Temperature (C)',           # Temperature (in Celsius)
                        'Pressure (Torr)',          # Pressure (in Torr)                        - Either low or high, since this can't really be controlled
                        'Flowrate (sccm)',          # Flowrate (in standard cubic centimeter)
                        'Precursor (mg)',           # Precursor mass (in milligrams)
                        'Precursor Ratio (:)']      # Ratio of two materials in the precursor 

    dependencies    = [ None,                       # The Temperature range depends on nothing
                        datanames[0],               # The pressure range depents on Temperature (i.e. Pressure is a function of P(T) )
                        None,                       # The Flowrate range depends on nothing
                        None,                       # The Precursor range depends on nothing
                        None ]                      # The Precursor ratio range depends on nothing.
    
    datalimits      = [ Temperature,                # An array must be used, as it is the independent variable for the Pressure
                        Pressure,                   # An array of the dependent function
                        [ 0, 200 ],                 # Min/Max of the flowrate
                        [ 100, 200 ],               # Min/Max of the Precursor amount
                        [ 0, 4 ] ]                  # Min/Max of the Precursor Ratio

    # Data types of each parameter.
    datatypes       = [ np.float32,
                        np.float32,
                        np.int8,
                        np.int8,
                        np.float32 ]

    sampler         = LHCS( samples       = samples, 
                          datanames       = datanames, 
                          dependencies    = dependencies, 
                          datalimits      = datalimits, 
                          datatypes       = datatypes, 
                          verbose         = True )

    sampler.run()
```
    This would create a sampling of data with a constraint on the pressure as a function of temperature
