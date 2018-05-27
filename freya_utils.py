import sys
import time

import numpy as np
import math as math
import json


################################################################################
# freya_read
################################################################################
def freya_read(file_name):

  file = open(file_name)
    
  initial_line = file.readline()
  initial_words = initial_line.split()

  if initial_words[1] is 253 or '253' or "253":
    initial_words[1] = 252
  
  lines = file.readlines()
  lines = lines[0:]
  file.close()
    
  import_begin = time.time()
    
  print('Importing Events...')

  nlines = len(lines)
  n = 1
  events = []
  event = []
  for i in range(0,nlines):
    line = lines[i]
    if( i == (nlines-1) ):
      events.append(event[:])
    elif ( not(line.split()[0].strip() == str(n+1)) ):
      event.append(line.strip())
    else:
      events.append(event[:])
      event = [line.strip()]
      n += 1
    
  # lines = [] #clear this portion of memory
    
  print(initial_words[3] + ' events imported.')
  import_end = time.time()
  import_time = import_end - import_begin
  print('Time: '+ str(import_time) + " sec")

  return events, initial_words

################################################################################
# freya_parse
################################################################################
def freya_parse(event_array, NEVENTSREAD, initial_words, En_cut):
  
    ranges_x = {} 

    ranges_x['A'] = [60,180]
    ranges_x['Z'] = [30,70]

    ranges_x['Product_A'] = ranges_x['A']
    ranges_x['Fragment_A'] = np.copy(ranges_x['Product_A'])

    ranges_x['n_Af'] = np.copy(ranges_x['Product_A'])
    ranges_x['m_Af'] = np.copy(ranges_x['Product_A'])
    ranges_x['TKE_A'] = np.copy(ranges_x['Product_A'])

    ranges_x['total_photon_energy'] = np.copy(ranges_x['Product_A'])
    ranges_x['energy_per_photon'] = np.copy(ranges_x['Product_A'])

    ranges_x['n_angular'] = [-1,1]

    ranges_x['n_mult'] = [0,20]
    ranges_x['m_mult'] = [0,40]

    ranges_x['n_TKE'] = [100,240]
    ranges_x['n_spectrum'] = [0,15]
    ranges_x['m_spectrum'] = [0,10]
    ranges_x['mannhart'] = [0.01, 100]
    ranges_x['n_A_TKE'] = [100, 240]


    ranges_y = {}

    ranges_y['Product_A'] = [0,12]
    ranges_y['Fragment_A'] = np.copy(ranges_y['Product_A'])

    ranges_y['n_Af'] = [0,5]
    ranges_y['m_Af'] = [0,8]
    ranges_y['TKE_A'] = [100,240]

    ranges_y['total_photon_energy'] = [ 0,12 ]
    ranges_y['energy_per_photon'] = [ 0,3.7 ]

    ranges_y['n_angular'] = [0.2,1.2]

    ranges_y['n_mult'] = [0,0.6]
    ranges_y['m_mult'] = [0,0.25]

    ranges_y['n_TKE'] = [0,12]
    ranges_y['n_spectrum'] = [0,0.55]
    ranges_y['m_spectrum'] = [0,1.3]
    ranges_y['mannhart'] = [0,0.7]
    ranges_y['n_A_TKE'] = np.copy(ranges_x['A'])

    ranges_z = {}

    ranges_z['n_A_TKE'] = [None, None]


    output_keys = ranges_x.keys()

    bin_number = {}

    bin_number['n_A_TKE'] = 70
    bin_number['n_TKE'] = 70
    bin_number['n_angular'] = 20
    bin_number['n_spectrum'] = 30
    bin_number['m_spectrum'] = 20



    bin_width = {}

    for key in bin_number.keys():
        bin_width[str(key)] = (ranges_x[str(key)][1] - ranges_x[str(key)][0]) / bin_number[str(key)]

    freya_output = {}
    
    Fragment_A = np.zeros((ranges_x['A'][1]- ranges_x['A'][0] + 1 , 3))
    Fragment_A[:,0] = range(ranges_x['A'][0],ranges_x['A'][1]+1)

    Product_A = np.copy(Fragment_A)

    Freya_Ah = []
    Freya_Al = []

    Freya_Z = np.zeros((ranges_x['Z'][1]- ranges_x['Z'][0] + 1 , 3))
    Freya_Z[:,0] = range(ranges_x['Z'][0] , ranges_x['Z'][1]+1)

    flightz = []

    fheavyz = []

    n_Af = np.copy(Fragment_A)

    m_Af = np.copy(Fragment_A)

    TKE_A = np.copy(Fragment_A)
    KE_A = []
    TKE = []
    KEl = []
    KEh = []

    TXE = []
    Ul = []
    Uh = []

    #  mannhart = np.zeros((len(mannhart_bins) , 3))
    #  mannhart[:,0] = mannhart_bins

    n_mult = np.zeros((ranges_x['n_mult'][1] - ranges_x['n_mult'][0] , 3))
    n_mult[:,0] = range(ranges_x['n_mult'][0],ranges_x['n_mult'][1])

    m_mult = np.zeros((ranges_x['m_mult'][1] - ranges_x['m_mult'][0] , 3))
    m_mult[:,0] = range(ranges_x['m_mult'][0],ranges_x['m_mult'][1] )

    energy_per_photon = np.zeros((ranges_x['A'][1]- ranges_x['A'][0]  , 4))
    energy_per_photon[:,0] = range(ranges_x['A'][0],ranges_x['A'][1] )

    total_photon_energy = np.copy(Fragment_A)

    #  the neutron and gamma spectrums, and the angular correlation will be calculated from the list which we will fill event by event. 
    neutrons = []
    light_neutrons = []
    heavy_neutrons = []
    n_spec = []
    light_n_spec = []
    heavy_n_spec = []

    n_ang = []
    n_ang_cut = []

    gammas = []
    light_gammas = []
    heavy_gammas = []
    m_spec = []
    light_m_spec = []
    heavy_m_spec = []

    #  the n_A_TKE bin columns need to be 'repeated' and 'tiled' so as to get every combinatorial pair
    n_A_TKE = np.zeros((bin_number['n_A_TKE'] * (ranges_x['A'][1] - ranges_x['A'][0]) , 4))
    n_A_TKE_energy_bins = np.arange(ranges_x['n_A_TKE'][0], ranges_x['n_A_TKE'][1], bin_width['n_A_TKE'] )
    n_A_TKE_mass_bins = np.arange(ranges_x['A'][0],ranges_x['A'][1] )
    n_A_TKE[:,0] = np.repeat(n_A_TKE_energy_bins,int( len(n_A_TKE) / len(n_A_TKE_energy_bins) ) )
    n_A_TKE[:,1] = np.tile( n_A_TKE_mass_bins , int( len(n_A_TKE)  / len(n_A_TKE_mass_bins ) ) )

    n_TKE = np.zeros((bin_number['n_TKE'],4))
    n_TKE[:,0] = np.arange(ranges_x['n_TKE'][0], ranges_x['n_TKE'][1], bin_width['n_TKE'])
    
    begin_parse = time.time()
    print('Parsing events...')
    for event in event_array[:NEVENTSREAD]:
        initial = []
        light = []
        heavy = []
        momenta_n = []
        momenta_m = []

        #  assign nn to be the number of neutrons for this event

        nn = int( event[0].split()[4] )

        

        if(nn > 0):
            nn_tf = 1
        else:
            nn_tf = 0
        initial = event[0:nn_tf+2]

        m = int( event[0].split()[5] )


        if (m>0):
            m_tf = 1
        else:
            m_tf = 0


        #  nnl is the integer number of neutrons coming from the light fragment
        nnl = int( event[nn_tf+2].split()[4] )

        #  ngl is the integer number of gammas coming from the light fragment
        ngl = int( event[nn_tf+2].split()[5] )
        #  if we have a positive number for nnl, set the nnl_tf to be 1 for true
        if(nnl > 0):
            nnl_tf = 1
        #  if we don't have any neutrons for the light fragment, set nnl_tf to 0 for false
        else:
            nnl_tf = 0
        if(ngl > 0):
            ngl_tf = 1
        else:
            ngl_tf = 0

        #  if we have (resp. don't have) neutrons for the fragments for an event we have (resp. don't have) additional lines to include in what we separate out as belonging to this event
        light = event[nn_tf+2:nn_tf+2+nnl_tf+ngl_tf+2]

        #  now we do the same thing for the heavy fragment
        if len(event) <= nn_tf+2+nnl_tf+ngl_tf+2:
            print(nn_tf+2+nnl_tf+ngl_tf+2)
            print(event)

        nnh = int( event[nn_tf+2+nnl_tf+ngl_tf+2].split()[4] )
        ngh = int( event[nn_tf+2+nnl_tf+ngl_tf+2].split()[5] )

        if(nnh > 0):
            nnh_tf = 1
        else:
            nnh_tf = 0
        if(ngh > 0):
            ngh_tf = 1
        else:
            ngh_tf = 0
        heavy = event[nn_tf+2+nnl_tf+ngl_tf+2:nn_tf+2+nnl_tf+ngl_tf+2+nnh_tf+ngh_tf+2]

        Afrag_l = int( light[0].split()[2] )
        Afrag_h = int( heavy[0].split()[2] )

        Freya_Ah.append(Afrag_h)
        Freya_Al.append(Afrag_l)

        Zfrag_l = int( light[0].split()[1] )
        Zfrag_h = int( heavy[0].split()[1] )

        #  write the appropriate values into the appropriate arrays
        #  the contents of each line for each events is detailed in the document: event_record.pdf

        Fragment_A[Afrag_l - ranges_x['A'][0] ][1] += 1
        Fragment_A[Afrag_h - ranges_x['A'][0] ][1] += 1

        n_Af[Afrag_l - ranges_x['A'][0] ][1] += nnl
        n_Af[Afrag_l - ranges_x['A'][0] ][2] += nnl**2

        n_Af[Afrag_h - ranges_x['A'][0] ][1] += nnh
        n_Af[Afrag_h - ranges_x['A'][0] ][2] += nnh**2

        m_Af[Afrag_l - ranges_x['A'][0] ][1] += ngl
        m_Af[Afrag_l - ranges_x['A'][0] ][2] += ngl**2

        m_Af[Afrag_h - ranges_x['A'][0] ][1] += ngh
        m_Af[Afrag_h - ranges_x['A'][0] ][2] += ngh**2

        Aprod_l = Afrag_l - nnl
        Aprod_h = Afrag_h - nnh

        neutrons.append(nn)
        neutrons.append(nnl)
        light_neutrons.append(nnl)
        neutrons.append(nnh)
        heavy_neutrons.append(nnh)

        gammas.append(m)
        gammas.append(ngl)
        light_gammas.append(ngl)
        gammas.append(ngh)
        heavy_gammas.append(ngh)
        
        Product_A[Aprod_l - ranges_x['A'][0]][1] += 1
        Product_A[Aprod_h - ranges_x['A'][0]][1] += 1

        Freya_Z[Zfrag_l - ranges_x['Z'][0]][1] += 1
        Freya_Z[Zfrag_h - ranges_x['Z'][0]][1] += 1

        fheavyz.append(Zfrag_h)
        flightz.append(Zfrag_l)

        ntot = nnl + nnh + nn
        n_mult[ntot][1] += 1

        mtot = ngl + ngh + m
        m_mult[mtot][1] += 1

        event_TKE = float( light[1].split()[0] ) + float( heavy[1].split()[0] )
        event_hKE = float( heavy[1].split()[0] )
        event_lKE = float( light[1].split()[0] )

        event_Ul = float( light[0].split()[3] )
        event_Uh = float( heavy[0].split()[3] )
        event_TXE = event_Ul + event_Uh

        TKE_A[Afrag_l - ranges_x['A'][0]][1] += event_TKE
        TKE_A[Afrag_l - ranges_x['A'][0]][2] += event_TKE**2

        TKE_A[Afrag_h - ranges_x['A'][0]][1] += event_TKE
        TKE_A[Afrag_h - ranges_x['A'][0]][2] += event_TKE**2

        KE_A.append([Afrag_l,event_lKE])
        KE_A.append([Afrag_h,event_hKE])

        KEl.append(event_lKE)
        KEh.append(event_hKE)
        TKE.append(event_TKE)

        TXE.append(event_TXE)
        Ul.append(event_Ul)
        Uh.append(event_Uh)

        energy_bin =  np.searchsorted(n_A_TKE[:,0],event_TKE)
        energy_bin = Afrag_l - ranges_x['A'][0] + energy_bin
        n_A_TKE[energy_bin,2] += nnl
        n_A_TKE[energy_bin,3] += nnl**2
        energy_bin = Afrag_h - ranges_x['A'][0] + energy_bin + 1
        n_A_TKE[energy_bin,2] += nnh
        n_A_TKE[energy_bin,3] += nnh**2
        energy_bin = int(98 - ranges_x['A'][0] + energy_bin + 1)
        n_A_TKE[energy_bin,2] += nn
        n_A_TKE[energy_bin,3] += nn**2


        n_TKE_bin = np.searchsorted(n_TKE[:,0], event_TKE)
        n_TKE[n_TKE_bin,1] += ntot
        n_TKE[n_TKE_bin,2] += ntot ** 2
        n_TKE[n_TKE_bin,3] += 1

        #  treat each case of neutrons and gamma separately (as determined by the 1 or 0 value of the indicating boolean nnl_tf for ex.)
        #  for each of these cases we write the appropriate values into the appropriate lists and arrays

        if(nn_tf > 0):
            parts = initial[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                light_n_spec.append(En)

        if( (nnl_tf > 0) and (ngl > 0) ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                light_n_spec.append(En)
            parts = light[3].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)
                light_m_spec.append(Eg)

                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][1] += Eg
                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][2] += Eg**2
                energy_per_photon[Afrag_l - ranges_x['A'][0] + 1][3] += 1

                total_photon_energy[Afrag_l - ranges_x['A'][0] + 1][1] += Eg
                total_photon_energy[Afrag_l - ranges_x['A'][0] + 1][2] += Eg**2
        elif( nnl_tf > 0 ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                light_n_spec.append(En)
        elif( ngl > 0 ):
            parts = light[2].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)
                light_m_spec.append(Eg)

                energy_per_photon[Afrag_l - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_l - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_l - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_l - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_l - ranges_x['A'][0] ][2] += Eg**2

        if( (nnh_tf > 0) and (ngh > 0) ):
            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                heavy_n_spec.append(En)
            parts = heavy[3].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)
                heavy_m_spec.append(Eg)

                energy_per_photon[Afrag_h - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_h - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
        elif( nnh_tf > 0 ):
            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                En = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_n.append( momentum )
                n_spec.append(En)
                heavy_n_spec.append(En)
        elif( ngh > 0 ):

            parts = heavy[2].split()
            for j in range( 0,int(len(parts)/4) ):
                Eg = float(parts[j*4])
                momentum = []
                momentum.append( float(parts[j*4 +1]) )
                momentum.append( float(parts[j*4 +2]) )
                momentum.append( float(parts[j*4 +3]) )
                momenta_m.append( momentum )
                m_spec.append(Eg)
                heavy_m_spec.append(Eg)

                energy_per_photon[Afrag_h - ranges_x['A'][0] ][1] += Eg
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][2] += Eg**2
                energy_per_photon[Afrag_h - ranges_x['A'][0] ][3] += 1

                total_photon_energy[Afrag_h - ranges_x['A'][0] ][1] += Eg
                total_photon_energy[Afrag_h - ranges_x['A'][0] ][2] += Eg**2

        for p in range(0,len(momenta_n)-1):
            for q in range(p+1,len(momenta_n)):
                dotprod = momenta_n[p][0] * momenta_n[q][0] + momenta_n[p][1] * momenta_n[q][1] + momenta_n[p][2] * momenta_n[q][2]
                n_ang.append(dotprod)
                if En >= En_cut: 
                    n_ang_cut.append(dotprod)

    begin_calc = time.time()

    with np.errstate(divide='ignore', invalid='ignore'):
        long_frag_counts = np.tile(Fragment_A[:,1][:-1],  int(len(n_A_TKE)/(ranges_x['A'][1] - ranges_x['A'][0])))
        n_A_TKE[:,2] = np.divide(n_A_TKE[:,2], long_frag_counts)
        n_A_TKE[:,3] = np.divide(n_A_TKE[:,3], long_frag_counts)
        n_A_TKE[:,3] = np.subtract(n_A_TKE[:,3],np.square(n_A_TKE[:,2]))
        n_A_TKE[:,3] = np.sqrt(n_A_TKE[:,3])
        n_A_TKE[:,2] = np.array( n_A_TKE[:,2], dtype = np.float)
        n_A_TKE[n_A_TKE == 0] = 'nan'

        total_photon_energy[:,1] = np.divide(total_photon_energy[:,1], Fragment_A[:,1])
        total_photon_energy[:,2] = np.divide(total_photon_energy[:,2], Fragment_A[:,1])
        total_photon_energy_squared = np.divide(total_photon_energy[:,1], Fragment_A[:,1])
        total_photon_energy[:,2] = np.subtract(total_photon_energy[:,2], total_photon_energy_squared)
        total_photon_energy[:,2] = np.sqrt(total_photon_energy[:,2])

        energy_per_photon[:,1] = np.divide(energy_per_photon[:,1], energy_per_photon[:,3])
        energy_per_photon[:,2] = np.divide(energy_per_photon[:,2], energy_per_photon[:,3])
        energy_per_photon[:,2] = np.subtract(energy_per_photon[:,2], np.square(energy_per_photon[:,1]))
        energy_per_photon[:,2] = np.sqrt(energy_per_photon[:,2])

        n_Af[:,1] = np.divide(n_Af[:,1] , Fragment_A[:,1])
        n_Af[:,2] = np.divide(n_Af[:,2] , Fragment_A[:,1])
        n_Af[:,2] = np.subtract(n_Af[:,2],np.square(n_Af[:,1]))
        n_Af[:,2] = np.sqrt(n_Af[:,2])

        TKE_A[:,1] = np.divide(TKE_A[:,1] , Fragment_A[:,1])
        TKE_A[:,2] = np.divide(TKE_A[:,2] , Fragment_A[:,1])
        TKE_A[:,2] = np.subtract(TKE_A[:,2],np.square(TKE_A[:,1]))
        TKE_A[:,2] = np.sqrt(TKE_A[:,2])
        

        m_Af[:,1] = np.divide(m_Af[:,1],Fragment_A[:,1])
        m_Af[:,2] = np.divide(m_Af[:,2],Fragment_A[:,1])
        m_Af[:,2] = np.subtract(m_Af[:,2],np.square(m_Af[:,1]))
        m_Af[:,2] = np.sqrt(m_Af[:,2])

        #  hist, bins = np.histogram(n_spec,bins = mannhart_bins)

        #  mannhart[:-1,1] = hist
        #  mannhart[:-1,2] = np.square(hist)
        #  mannhart[:,1] = np.divide(mannhart[:,1],len(n_spec))
        #  mannhart[:,1] = np.divide(mannhart[:,1],mannhart_bindiff.flatten())
        #  mannhart[:,2] = np.divide(mannhart[:,2],sum(n_spec))
        #  mannhart[:,2] = np.subtract(mannhart[:,2],np.square(mannhart[:,1]))
        #  mannhart[:,2] = np.sqrt(mannhart[:,2])
        #  mannhart[:,2] = np.divide(mannhart[:,2],mannhart_bindiff.flatten() * 1000  )

        n_TKE[:,1] = np.divide(n_TKE[:,1], n_TKE[:,3])
        n_TKE[:,2] = np.divide(n_TKE[:,2], n_TKE[:,3])
        n_TKE[:,2] = np.subtract(n_TKE[:,2], np.square(n_TKE[:,1]))
        n_TKE[:,2] = np.sqrt(n_TKE[:,2])

        nubar = float(len(n_spec))/float(len(event_array))
        mubar = float(len(m_spec))/float(len(event_array))

        nu1 = nubar 
        nu2 = 0.0
        nu3 = 0.0
        nu4 = 0.0

        n_mult[:,1] = np.divide(n_mult[:,1] , np.sum(n_mult[:,1]) )

        m_mult[:,1] = np.divide(m_mult[:,1] , np.sum(m_mult[:,1]) )

        for i in n_mult[:,0]:
            i = int(i)
            if sum(n_mult[:,1]) != 0:
                nu2 += (i) * (i-1) * n_mult[i,1]
                nu3 += (i) * (i-1) * (i-2) *  n_mult[i,1]
                nu4 += (i) * (i-1) * (i-2) * (i-3) * n_mult[i,1]

        end_calc = time.time()

    print('Time: ' + str(end_calc - begin_calc) + ' sec')

    number_of_bins = int(np.ceil(bin_number['n_angular']))
    hist, bin_edges = np.histogram(n_ang, bins = number_of_bins, range = (ranges_x['n_angular'][0],ranges_x['n_angular'][1]), normed = True)
    n_angular = np.zeros((number_of_bins,3))
    n_angular[0] = [ -1 , 1 , bin_width['n_angular']]
    n_angular[:,0] = np.add(bin_edges, bin_width['n_angular']/2)[:-1]
    n_angular[:,1] = hist
    n_angular[:,2] = np.nan

    number_of_bins = int(np.ceil(bin_number['n_angular']))
    hist, bin_edges = np.histogram(n_ang_cut, bins = number_of_bins, range = (ranges_x['n_angular'][0],ranges_x['n_angular'][1]), normed = True)
    n_angular_cut = np.zeros((number_of_bins,3))
    n_angular_cut[0] = [ -1 , 1 , bin_width['n_angular']]
    n_angular_cut[:,0] = np.add(bin_edges, bin_width['n_angular']/2)[:-1]
    n_angular_cut[:,1] = hist
    n_angular_cut[:,2] = np.nan

    n_spectrum_bin_widths = bin_width['n_spectrum']
    number_of_bins = int(np.ceil(bin_number['n_spectrum']))
    hist, bin_edges = np.histogram(n_spec, bins = number_of_bins, range = (ranges_x['n_spectrum'][0],ranges_x['n_spectrum'][1]), normed = True)
    n_spectrum = np.zeros((number_of_bins ,3))
    n_spectrum[0] = [0,ranges_x['n_spectrum'][1],n_spectrum_bin_widths]
    n_spectrum[:,0] = np.add(bin_edges,n_spectrum_bin_widths/2)[:-1]
    n_spectrum[:,1] = hist
    n_spectrum[:,2] = np.nan

    m_spectrum_bin_widths = bin_width['m_spectrum']
    number_of_bins = int(np.ceil(bin_number['m_spectrum']))
    hist, bin_edges = np.histogram(m_spec, bins = number_of_bins, range = (ranges_x['m_spectrum'][0] - m_spectrum_bin_widths, ranges_x['m_spectrum'][1]), normed = True)
    m_spectrum = np.zeros((number_of_bins,3))
    m_spectrum[0] = [0,ranges_x['m_spectrum'][1],m_spectrum_bin_widths]
    m_spectrum[:,0] = np.add(bin_edges,m_spectrum_bin_widths/2)[:-1]
    m_spectrum[:,1] = hist
    m_spectrum[:,2] = np.nan

    Fragment_A[:,1] = np.divide(Fragment_A[:,1], 2 * NEVENTSREAD)

    Product_A[:,1] = np.divide(Product_A[:,1], 2 * NEVENTSREAD)

    Freya_Z[:,1] = np.divide(Freya_Z[:,1], 2 * NEVENTSREAD)


    freya_output['Z00'] = initial_words[0]
    freya_output['A00'] = initial_words[1]
    freya_output['isotope'] = 'Z = ' + str(freya_output['Z00']) + ' A = ' + str(freya_output['A00'])
    freya_output['number_of_events'] = initial_words[3]
    freya_output['nubar'] = nubar
    freya_output['mubar'] = mubar
    freya_output['nu1'] = nu1
    freya_output['nu2'] = nu2
    freya_output['nu3'] = nu3
    freya_output['nu4'] = nu4
    freya_output['Fragment_A'] = Fragment_A
    freya_output['Product_A'] = Product_A
    freya_output['n_A_TKE'] = n_A_TKE
    #  freya_output['mannhart'] = mannhart
    freya_output['n_angular'] = n_angular
    freya_output['n_angular_cut'] = n_angular_cut
    freya_output['n_ang'] = n_ang

    freya_output['neutrons'] = np.array(neutrons)
    freya_output['light_neutrons'] = np.array(light_neutrons)
    freya_output['heavy_neutrons'] = np.array(heavy_neutrons)
    freya_output['n_mult'] = n_mult

    freya_output['gammas'] = np.array(gammas)
    freya_output['light_gammas'] = np.array(light_gammas)
    freya_output['heavy_gammas'] = np.array(heavy_gammas)
    freya_output['all_gammas'] = np.add(freya_output['light_gammas'], freya_output['heavy_gammas'])
    freya_output['m_mult'] = m_mult

    freya_output['n_spec'] = np.array(n_spec)
    freya_output['light_n_spec'] = np.array(light_n_spec)
    freya_output['heavy_n_spec'] = np.array(heavy_n_spec)
    freya_output['n_spectrum'] = n_spectrum

    freya_output['m_spec'] = np.array(m_spec)
    freya_output['light_m_spec'] = np.array(light_m_spec)
    freya_output['heavy_m_spec'] = np.array(heavy_m_spec)
    freya_output['m_spectrum'] = m_spectrum

    freya_output['total_photon_energy'] = total_photon_energy
    freya_output['energy_per_photon'] = energy_per_photon
    freya_output['n_Af'] = n_Af
    freya_output['TKE_A'] = TKE_A
    freya_output['KE_A'] = np.array(KE_A)
    freya_output['TKE'] = np.array(TKE)
    freya_output['KEl'] = np.array(KEl)
    freya_output['KEh'] = np.array(KEh)
    freya_output['m_Af'] = m_Af
    freya_output['n_TKE'] = n_TKE
    freya_output['Ah'] = Freya_Ah
    freya_output['Al'] = Freya_Al
    freya_output['Z'] = Freya_Z
    freya_output['Z_l'] = np.array(flightz)
    freya_output['Z_h'] = np.array(fheavyz)
    freya_output['Ul'] = (np.array(Ul))
    freya_output['Uh'] = (np.array(Uh))
    freya_output['TXE'] = (np.array(TXE))

    print('All events parsed.')
    end_parse = time.time()
    parse_time = end_parse - begin_parse
    print('Time: '+ str(parse_time) + ' sec')

    return freya_output, ranges_x

