'''
    Name:           Latin Hypercube Sampler class for dynamic data

    Description:    Takes an input of the CVD setup and provides initial guesses to maximimize stastical difference between tests

    Date:           18 October 2022

    Author:         John Ferrier, NEU Physics
'''

from math import ceil
import os
import random
import numpy as np
import pandas as pd
from scipy.stats import _qmc
import matplotlib.pyplot as plt
from Interpolator import Interpolator


class LHCS:

    # Initialization function
    def __init__( self,
                    samples             = 20,               # The amount of LHC samples wanted
                    datanames           = [],               # The names of the data categories
                    dependencies        = [],               # What each data option depends on
                    datalimits          = [],               # The min and max value range of the data (min and max can be an array of the same length)
                    datatypes           = [],               # The data type to be used (int or float)
                    verbose             = False  ) -> None: # Whether or not to print live information or not
        
        self.samples        = samples
        self.datanames      = datanames
        self.datalimits     = datalimits
        self.dependencies   = dependencies
        self.datatypes      = datatypes
        self.verbose        = verbose

        self.breakInitial   = False                         # Whether or not the stop the initialization process

        # Check that none of the values are empty and that all values match in length
        data_length         = len( datanames )
        if not len( self.dependencies ) == data_length:
            if self.verbose:
                print( "Dependency list length does not match data provided!\nSetting all to `None`" )
            # Just set this as full None, since no dependencies doesn't mean this can't work
            self.dependencies   = [ None for i in range( data_length ) ]
        
        if not len( self.datatypes ) == data_length:
            if self.verbose:
                print( "Datatypes list length does not match data provided!\nSetting all to `float64`" )
            # If not set, just presume float
            self.datatypes  = [ np.float64 for i in range( data_length ) ]

        if not len( self.datalimits ) == data_length:
            print( "Datalimits list length does not match data provided!\nStopping Sampler" )
            # This legitimately needs to match
            self.breakInitial   = True
        
        if not self.breakInitial:

            # Check for dependencies
            self.isDependent    = False                 # Whether there are dependencies or not
            self.Dependent_idx  = []                    # The indice combos
            self.dependent_lst  = []                    # All indices used (in case of multiples)

            if self.verbose:
                print( "Checking Dependencies..." )

            for i, d in enumerate( self.dependencies ):
                if d is not None:
                    # Find the index of the dependent and independent values
                    if d in self.datanames:
                        self.Dependent_idx.append( [ i, self.datanames.index(d) ] )
                        print( f"`{self.datanames[i]}` depends on `{d}`" )
                        self.dependent_lst.append( i )
                        self.dependent_lst.append( self.datanames.index(d) )
                    else:
                        self.breakInitial = True
                        print( "The independent value does not exist!\nPlease check your spelling!" )

            # Remove redundancies
            if len( self.dependent_lst ) > 0:
                self.dependent_lst  = [ *set( self.dependent_lst ) ]

            # Ensure that each dependency has the same amount of data points that correspond.
            if len( self.Dependent_idx ) > 0 and not self.breakInitial:
                if self.verbose:
                    print( "Comparing dependency lengths..." )

                for d in self.Dependent_idx:

                    if not len( self.datalimits[d[0]] ) == len( self.datalimits[d[1]] ):
                        print( f"`{self.datanames[d[0]]}` data length does not match the data length for `{self.datanames[d[1]]}`!\nStopping sampler" )
                        self.breakInitial   = True
    
    # Makes guess for experiments based on inputs
    def run( self ):

        # Check the ensure that the sampler can be run
        if not self.breakInitial:

            # Make initial guess for the non-dependent items
            self.lhc        = _qmc.LatinHypercube(  d               = len( self.datanames ) - len( self.dependent_lst ),    # Dimensions of the hypercube minus the dependency values
                                                    centered        = True,                                                 # Usually better for CVD setups
                                                    optimization    = "random-cd" )                                         # Increased computational cost but it lowers the discrepancy value

            self.guess     = self.lhc.random( n = self.samples )
            
            # Make initial guess for the dependent items
            self.man_guess = self.manualGuess( showPlot = True )

            # Scales the outputs to the min/max values
            self.scaleOutputs()

            # Set type
            #self.setDataType()

            self.guess  = self.guess.tolist()

            self.dependent_lst.sort()

            for i in range( len( self.guess ) ):

                for j, mg in enumerate( self.man_guess[i] ):

                    self.guess[i].insert( self.dependent_lst[j], self.man_guess[i][j] )

            # Output the data
            self.df = pd.DataFrame( self.guess, columns = self.datanames )
            print( self.df.to_markdown() )

    # Scales the non-dependent values
    def scaleOutputs( self ):

        # Remove dependent parts
        self.dependent_lst.sort( reverse = True )       # Reverse the order to avoid issues with removing points
        for d in self.dependent_lst:

            del self.datalimits[ d ]
            del self.datatypes[ d ]

        # Get min/max values of non-dependent values

        minimum         = np.array( [ m[0] for m in self.datalimits ] )
        maximum         = np.array( [ m[1] for m in self.datalimits ] )

        # Scale outputs
        self.guess      = _qmc.scale( self.guess, l_bounds = minimum, u_bounds = maximum )

    # Forces the dataype on the non-dependent values
    def setDataType( self ):

        #### TOO MANY ERRORS FROM THIS. SKIP FOR NOW

        # Cycle through each dimension
        for i,t in enumerate( self.datatypes ):

            tmp                     = self.guess[ :, i ].copy()
            self.guess[ :, i ]      = tmp.astype( t )

     
    # Make a manual guess
    def manualGuess( self, showPlot = False):

        if self.verbose:
            print( "Starting the manual guess..." )

        # Cycle through each dependency
        for d in self.Dependent_idx:

            # Samples to consider
            N       = self.samples

            init_x  = self.datalimits[ d[1] ]
            init_y  = self.datalimits[ d[0] ]

            # Range Resolution
            x_min   = init_x[ 0 ]
            x_max   = init_x[ -1 ]
            y_min   = init_y[ 0 ]
            y_max   = init_y[ -1 ]

            Δx      = x_max - x_min
            Δy      = y_max - y_min

            # Get finite integral of f(x)
            tmp_stp = 5000                                          # Arbitrary resolution of a fit

            # If less data than 5000 points
            if len( init_x ) < tmp_stp:
                temp_x, temp_y  = self.interpolate( init_x, init_y, tmp_stp )
            else:
                tmp_stp     = len( init_x )
                temp_x      = init_x
                temp_y      = init_y

            δx      = Δx/tmp_stp
            

            stp_mul = 0.
            N_stps  = 0

            # Start low and increase to reduce computational time. This loop is easy to do
            while N_stps < 2:
                
                stp_mul += 1.

                stps    = stp_mul*N               # steps chosen for sampling purposes. Will increase to insure enough variance.

                # Amount of sample squares (under the curve)
                N_s     = stps*stps

                # Area under the curve
                int_y   = δx*np.sum( temp_y )

                # Area of rectangular space
                A       = Δx*Δy

                # Ratio of undercurve vs entire rectangle
                N_r     = int_y/A 

                # Amount of squares needed for entire ΔxΔy rectangle
                Squares     = ceil( N_s/N_r )

                # scaling constant for the squares such that Squares = C^2 * A
                C           = np.sqrt( Squares/A )
                N_x         = round( C*Δx )
                N_y         = round( C*Δy )

                # See if ratios follow Total should be more or equal to N_s
                Total   = round( N_x*N_y*N_r )

                # Created centered mesh grid
                x_shift     = Δx/( 2*N_x )
                x           = np.linspace( x_shift, x_max - x_shift, N_x )

                y_shift     = Δy/( 2*N_y )
                y           = np.linspace( y_shift, y_max - y_shift, N_y )

                # Ratio of needed squares (N_s) versus the amount that actually exist under the curve
                N_stps  = float(Total)/(N**2)

                print( N_stps )
            #j = 1

        
            ### Steps are greater than 2. Move on  ###


            # Cycle through each point
            options_list = []
            for i, xval in enumerate( x ):

                # Get max y (from the function) value for this x value
                # Find closest x-value in temp_x
                x_idx       = np.abs( temp_x - xval ).argmin()
                yval_max = temp_y[x_idx]

                # Get the y value max indice
                idx     = ( np.abs(y - yval_max) ).argmin()

                # Cycle through these indices, setting them equal to 1
                for j in range( idx + 1 ):
                    options_list.append( [i, j] )

            # Choose random points
            selected_pnts   = []

            for k in range( N ):

                #print( f"Choices available = {len(options_list)}" )
                choice  = random.choice( options_list )
                #print( choice )
                selected_pnts.append( choice )

                # Remove those indices
                options_list = self._remove_indices( options_list, choice[0], choice[1] )


            ret     = [ [ x[ p[0] ], y[ p[1] ] ] for p in selected_pnts ]

            if showPlot:
                # Build points
                x_pnts  = [ x[p[0]] for p in selected_pnts ]
                y_pnts  = [ y[p[1]] for p in selected_pnts ]

                # Now, plot the points
                fig = plt.figure()
                ax  = fig.add_subplot(1,1,1)

                minorx_ticks     = np.linspace( 0, x_max, N_x+1 )
                minory_ticks     = np.linspace( 0, y_max, N_y+1 )

                ax.set_xticks( minorx_ticks, minor = True )
                ax.set_yticks( minory_ticks, minor = True )

                ax.grid( which = 'both' )
                ax.grid( which = 'minor', alpha = 0.5 )

                ax.plot( temp_x, temp_y )
                ax.plot( x_pnts, y_pnts, 'o', color = 'red' )
                plt.show()
        
        return ret

    # Removes indices from lists when handling manual dependencies
    def _remove_indices( self, lst, i, j ):
    
        matching = []

        for k in range( len( lst ) ):

            if lst[k][0] == i:
                matching.append( k )
            if lst[k][1] == j:
                matching.append( k )

        # remove redundant indices
        matching = [*set(matching)]

        # Sort in descending order
        matching.sort( reverse = True )

        # Remove values from list
        for l in matching:
            del lst[l]

        # Return remainder
        return lst

    # Fits a function given a dependency
    def interpolate( self, x, y, points = 1000, showplot = False ):
        
        self.INTER  = Interpolator( x, y, points )

        if showplot:
            self.INTER.plot_Old_v_New()

        return self.INTER.new_x, self.INTER.new_y


if __name__ == "__main__":

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

    sampler     = LHCS( samples         = samples, 
                        datanames       = datanames, 
                        dependencies    = dependencies, 
                        datalimits      = datalimits, 
                        datatypes       = datatypes, 
                        verbose         = True )

    #sampler.interpolate( x = Temperature, y = Pressure, showplot = True )

    sampler.run()
    #sampler.save( os.path.join( os.getcwd(), f'LHCS_IntialSamples_x{samples}.csv' ) )