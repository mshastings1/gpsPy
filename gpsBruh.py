#!/usr/bin/env python3.7

class gpsBruh:
    '''
    author: Mitchell Hastings
    Origin Date: 07/16/2021
    Class Description:
        > This is a class to deal with gps data from University of Nevada Reno
        > MAGNET program. It's going to host methods for reading in the data
        > from the complicated file formats riddled with extra delimiters
        > and be able to store the time series, station location, model parameters
        > rms errors, and whatever else I damn well wish
    Class Attributes:
        > data                  - time series dataframe of all componenets (fills out with some methods)
        > site                  - station name
        > lon, lat, elev        - station longitude, station latitude, station elevation
        > eqDates               - yymmmdd of the earthquakes in the series
        > eqTimes               - decimal year of the eqDates
        > model                 - dictionary of model estimation from weighted least squares
            > Contains 3 subdictionaries for North, East, and Up components.
            > e.g. use model.get('north') to get the model parameters dictionary for north
            > e.g. use gpsBruh.rms() to compute rms for all models and add the value to dictionary
        > M_N, M_E, M_U         - matrix forms of model estimations, useful for computing forward models
        > more shit to add when my virtual overlords request it
    Class Methods:
        > methods()                  - list all available methods for the class and their descriptions
        > load_UNR(file)             - Read in data from the UNR file types
        > time_window(start, stop)   - subset the data between two decimal year times 
        > read_eqDates(file)         - read a file with yymmmdd dates for eqs and add them to the object
        > get_eqTimes()              - convert eqDates to decimal years
        > neu_wls()                  - weighted least squares for north, east, and up components
        > forward_model(t,m)         - forward model of my fancy algorithm for a vectorized time array
                                       and model parameter set. 
        > model_gps()                - run my fabulous weighted least squares algorithm and compute forward 
                                       models for each component to be used in plotting
        > save_pickleBruh()          - save the object to a pickle file 
        > plotNEU()                  - plot the North, East, Up componenets
    Supporting Functions:
        > load_pickleBruh(file)      - load a pickled file of a gpsBRUH instance. the __init__ doesn't like 
                                       trying to load a pickled file because of conflicts with local variables.
    '''

    def __init__(self, *args):
        
        #############################################
        ##                                         ##
        ##         The Almighty Constructor!       ##
        ##                                         ##
        #############################################

        print("#--- WELCOMETOTHEHOGCRANKERSBALL ---#")

        if (len(args) > 1):
            raise ValueError("TOO MANY ARGS, COWBOY. JUST GIVE FILENAME OR NO ARGS")
        elif (len(args) == 1):
            if (args[0][-4:] == '.pkl'):
                print("ILLEGAL INITIALIZATION, PARTNER! use \'load_pickleBRUH\' function to load a pickle file")
            else:
                self.load_UNR(args[0])
        else:
            print("# CREATING EMPTY GPSBRUH OBJ #")
        
    def load_UNR(self, filename):

        from pandas import DataFrame
        from os import path

        # check to see if file exists and read it in if it does exist
        if (path.exists(filename) == False):
            raise ValueError("I CANNOT FIND %s PARTNER" % (filename))
        else:
            print("READING FILE, COWBOY: %s" % (filename))
            f = open(filename,'r')

        # initialize lists for all the columns, probs a better way to do this with positional lists and comprehensions
        s_l = []; y_l = []; dy_l = []; m_l = []; w_l = []; d_l = []; r_l = [];
        e0_l = []; e_l = []; n0_l = []; n_l = []; u0_l = []; u_l = []; a_l = [];
        se_l = []; sn_l = []; su_l = []; cen_l = []; ceu_l = []; cnu_l = [];

        # loop over the file
        while True:

            # read the first line (and then subsequent lines through the loop)
            line = f.readline()

            # exit when end of file is found
            if line == '':
                break

            # break up line on spaces
            split = line.split(' ')

            # if it is the header line - skip it
            if split[0] == 'site':
                continue

            # loop over the split string of data and remove empty spaces for pushing
            i = 0                           # index variable
            list_end = len(split)           # initialize how many items in list
            while (i != list_end):          # loop over each index and remove empty elements
                if (split[i] == ''):
                    split.pop(i)
                    list_end -= 1
                    i -= 1
                else:
                    i += 1

            # append elements to list 
            s_l.append(split[0]); y_l.append(split[1]); dy_l.append(float(split[2])); m_l.append(int(split[3]));
            w_l.append(int(split[4])); d_l.append(int(split[5])); r_l.append(float(split[6])); e0_l.append(float(split[7]));
            e_l.append(float(split[8])); n0_l.append(float(split[9])); n_l.append(float(split[10])); u0_l.append(float(split[11]));
            u_l.append(float(split[12])); a_l.append(float(split[13])); se_l.append(float(split[14])); sn_l.append(float(split[15]));
            su_l.append(float(split[16])); cen_l.append(float(split[17])); ceu_l.append(float(split[18])); cnu_l.append(float(split[19]));

        # create a dataframe object to receive data
        df = DataFrame(data=None, columns=['YYMMMDD','decYear','MJD','week','day','reflon','e0','east','n0','north','u0','up','ant',\
                'sig_e','sig_n','sig_u','corr_en','corr_eu','corr_nu'])

        # set lists to the dataframe
        df.YYMMMDD = y_l; df.decYear = dy_l; df.MJD = m_l; df.week = w_l; df.day = d_l; df.reflon = r_l; df.e0 = e0_l; df.east = e_l;
        df.n0 = n0_l; df.north = n_l; df.u0 = u0_l; df.up = u_l; df.ant = a_l; df.sig_e = se_l; df.sig_n = sn_l; df.sig_u = su_l; df.corr_en = cen_l;
        df.corr_eu = ceu_l; df.corr_nu = cnu_l;

        self.data = df
        self.site = s_l[0]

    def time_window(self, start, stop):
        
        newframe = self.data[(self.data.decYear >= start) & (self.data.decYear <= stop)]
        newframe = newframe.reset_index(); newframe = newframe.drop(axis=1,labels='index')
        self.data = newframe

        print("#--- FLOURIDEINTHEHARDSELTZERWATER ---#")


    def plotNEU(self):

        from matplotlib.pyplot import figure

        fig = figure()
        axU = fig.add_subplot(313)
        axN = fig.add_subplot(311)
        axE = fig.add_subplot(312)

        axN.plot(self.data.decYear, (self.data.north*1e3)-(self.data.north[0]*1e3), 'bo', ms=1)
        axE.plot(self.data.decYear, (self.data.east*1e3)-(self.data.east[0]*1e3), 'bo', ms=1)
        axU.plot(self.data.decYear, (self.data.up*1e3)-(self.data.up[0]*1e3), 'bo', ms=1)

        axN.set_xticks([])
        axE.set_xticks([])
        
        axU.set_xlabel('Year')
        axN.set_ylabel('North (mm)')
        axE.set_ylabel('East (mm)')
        axU.set_ylabel('Up (mm)')
        axN.set_title(self.site)

        print("#--- MAYBANGERSRESIDEINYOURHEART ---#")

        return fig

    def read_eqDates(self, filename):

        f = open(filename,'r')

        eqs = []

        while True:
            line = f.readline()
            if line == '':
                break
            date = line.strip('\n')
            eqs.append(date)

        self.eqDates = eqs

        print("# LIFTHEAVYSTONEMAKESADHEADVOICEQUIET #")

    def get_eqTimes(self):

        num_eqs = len(self.eqDates)

        eq_in_series = []
        new_eqDates = []

        for i in range(0,len(self.data)):
            for e in self.eqDates:
                if (self.data.YYMMMDD[i] == e):
                    eq_in_series.append(self.data.decYear[i])
                    new_eqDates.append(self.data.YYMMMDD[i])

        self.eqTimes = eq_in_series
        self.eqDates = new_eqDates

        print("# MIDWESTLAWNCAREDADSWHOSMASHBREWS #")

    def neu_wls(self):
        '''
        weighted least squares solution for modeling annual and semiannual signals
        '''

        # import necessary functions and modules
        from numpy import sin, cos, pi, matmul, matlib, linalg, dot

        # number of eqs to add
        numEQ = len(self.eqTimes)

        # modeling frequencies (cycles per year)
        self.annual_freq = 1; self.semiannual_freq = 2

        # create the G matrix 
        G = matlib.empty((len(self.data.decYear),6+numEQ))
        for i in range(0,len(self.data.decYear)):
            G[i,0] = 1
            G[i,1] = self.data.decYear[i]
            G[i,2] = sin(2*pi*self.annual_freq*self.data.decYear[i])    # annual
            G[i,3] = cos(2*pi*self.annual_freq*self.data.decYear[i])    # annual
            G[i,4] = sin(2*pi*self.semiannual_freq*self.data.decYear[i])    # semi-annual
            G[i,5] = cos(2*pi*self.semiannual_freq*self.data.decYear[i])    # semi-annual
            for j in range(0,numEQ):
                if (self.data.decYear[i] <= self.eqTimes[j]):
                    G[i,6+j] = 0
                else:
                    G[i,6+j] = 1

        # construct the weight matrices
        Wn = matlib.eye(len(self.data)); We = matlib.eye(len(self.data)); Wu = matlib.eye(len(self.data))
        for i in range(0,len(self.data)):
            Wn[i,i] = 1/self.data.sig_n[i]
            We[i,i] = 1/self.data.sig_e[i]
            Wu[i,i] = 1/self.data.sig_u[i]

        # disps relative to first point
        n_disp = self.data.north - self.data.north[0]
        e_disp = self.data.east - self.data.east[0]
        u_disp = self.data.up - self.data.up[0]

        # add disps to dataframe
        self.data['n_disp'] = n_disp
        self.data['e_disp'] = e_disp
        self.data['u_disp'] = u_disp

        # make lists to zip and loop through
        D = [n_disp, e_disp, u_disp]
        W = [Wn, We, Wu]

        # list to store parameter estimation
        M = []

        # first half of the least squares (up to multiplying by the data)
        G_t = G.transpose()
        for w,d in zip(W,D):
            a = matmul(G_t,w); b = matmul(a,G)
            c = linalg.inv(b); e = matmul(c,G_t)
            f = matmul(e,w); m = dot(f,d)
            M.append(m)

        # store matrices for use in forward model
        self.M_N = M[0]; self.M_E = M[1]; self.M_U = M[2]

        # creat dictionary to store model info in readable way
        model = {'north': {}, 'east': {}, 'up': {}}

        model['north']['y-int'] = M[0][0,0]
        model['north']['vel'] = M[0][0,1]
        model['north']['annual_sin'] = M[0][0,2]
        model['north']['annual_cos'] = M[0][0,3]
        model['north']['semiannual_sin'] = M[0][0,4]
        model['north']['semiannual_cos'] = M[0][0,5]
        
        model['east']['y-int'] = M[1][0,0]
        model['east']['vel'] = M[1][0,1]
        model['east']['annual_sin'] = M[1][0,2]
        model['east']['annual_cos'] = M[1][0,3]
        model['east']['semiannual_sin'] = M[1][0,4]
        model['east']['semiannual_cos'] = M[1][0,5]

        model['up']['y-int'] = M[2][0,0]
        model['up']['vel'] = M[2][0,1]
        model['up']['annual_sin'] = M[2][0,2]
        model['up']['annual_cos'] = M[2][0,3]
        model['up']['semiannual_sin'] = M[2][0,4]
        model['up']['semiannual_cos'] = M[2][0,5]

        c = 6
        for key in self.eqDates:
            model['north']['EQ_'+key] = M[0][0,c]
            model['east']['EQ_'+key] = M[1][0,c]
            model['up']['EQ_'+key] = M[2][0,c]
            c+=1

        self.model = model

        print("#--- POSTALMONDCLARITY ---#")


    def forward_model(self,t,m):

        from numpy import sin, cos, pi

        y = m[0,0] + m[0,1]*t + m[0,2]*sin(2*pi*t*self.annual_freq) + m[0,3]*cos(2*pi*t*self.annual_freq) + m[0,4]*sin(2*pi*t*self.semiannual_freq) + m[0,5]*cos(2*pi*t*self.semiannual_freq)

        for i in range(0,len(self.eqTimes)):
            if (t >= self.eqTimes[i]):
                y += m[0,6+i]

        return y

    def rms(self):
        '''
        root mean square error 
        '''

        n = (self.data.n_disp - self.data.forward_north)**2
        e = (self.data.e_disp - self.data.forward_east)**2
        u = (self.data.u_disp - self.data.forward_up)**2

        n = (n.sum()/len(self.data))**0.5
        e = (e.sum()/len(self.data))**0.5
        u = (u.sum()/len(self.data))**0.5

        self.model['north']['rmse'] = n
        self.model['east']['rmse'] = e
        self.model['up']['rmse'] = u

        print("# THEMACHOMANbilmuri #") 
    
    def model_gps(self):
        '''
        forward model the least squares inversion
        '''

        # run least squares
        self.neu_wls()

        # initialize lists to append to
        n = []; e = []; u = []

        for i in range(0,len(self.data.decYear)):
            y_n = self.forward_model(self.data.decYear[i],self.M_N)
            n.append(y_n)

            y_e = self.forward_model(self.data.decYear[i],self.M_E)
            e.append(y_e)

            y_u = self.forward_model(self.data.decYear[i],self.M_U)
            u.append(y_u)


        # push forwards back to input dataframe
        self.data['forward_north'] = n
        self.data['forward_east'] = e
        self.data['forward_up'] = u
        
        # compute rms
        self.rms()

        print("# GITSHRECKED #") 

    def plotModel(self):

        from matplotlib.pyplot import figure

        fig = figure()
        axU = fig.add_subplot(313)
        axN = fig.add_subplot(311)
        axE = fig.add_subplot(312)

        for eq in self.eqTimes:
            axN.axvline(eq)
            axE.axvline(eq)
            axU.axvline(eq)

        axN.plot(self.data.decYear, self.data.n_disp*1e3, 'bo', ms=1)
        axE.plot(self.data.decYear, self.data.e_disp*1e3, 'bo', ms=1)
        axU.plot(self.data.decYear, self.data.u_disp*1e3, 'bo', ms=1)

        axN.plot(self.data.decYear, self.data.forward_north*1e3, 'r-', lw=1.5)
        axE.plot(self.data.decYear, self.data.forward_east*1e3, 'r-', lw=1.5)
        axU.plot(self.data.decYear, self.data.forward_up*1e3, 'r-', lw=1.5)

        axN.set_xticks([])
        axE.set_xticks([])

        axU.set_xlabel('Year')
        axN.set_ylabel('North (mm)')
        axE.set_ylabel('East (mm)')
        axU.set_ylabel('Up (mm)')
        axN.set_title(self.site)

        print("# MYTOP10MOSTBROOTALBREAKDOWNSOF2047 #") 

        return fig

    def save_pickleBruh(self):

        import pickle
        import os

        filename = './' + self.site + '.pkl'

        c = 0
        while (os.path.exists(filename)):
            filename = './' + self.site + str(c) + '.pkl'
            c += 1
            
        with open(filename, 'wb') as file:

            pickle.dump(self, file)

        print("MACHO MAN bilmuri CREATED FILE: %s" % (filename))
            
def load_pickleBruh(filename):

    import pickle

    with open(filename,'rb') as file:
         new_obj = pickle.load(file)

    return new_obj
