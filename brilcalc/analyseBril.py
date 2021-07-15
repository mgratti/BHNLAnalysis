import csv
import os
import pandas as pd


# 1. 'lumiEach' :    Get the luminosity associated to each specified path for the entire dataset and divided in periods (A,B,C,D)
# 2. 'lumiOverlap':  Get the overlap luminosity between specified trigger paths
# 3. 'lumiTot':      Get the luminosity associated to the sum of the specified trigger paths, after removing overlaps
# 4. 'lumiTotByPart' samee as lumiTot, but done only for specified "part" (more trigger parts were active at the same time, datasets are divided in parts)



# Note: All luminosities are in units of 1/fb

if __name__ == "__main__":

  import argparse

  parser = argparse.ArgumentParser(description='lumi calculations for b-parking')
  parser.add_argument('--whichana', type=str, dest='whichana', help='whichana', default='lumiEach', choices=['lumiEach', 'lumiTot', 'lumiTotByPart', 'lumiOverlap'])
  parser.add_argument('--part', type=str, dest='part', help='Part of the dataset you want to investigate, e.g. part0,part1,...,part5', 
                                choices=['part0', 'part1', 'part2', 'part3', 'part4', 'part5' ], default=None)

  opt = parser.parse_args() 

  #### DO NOT EDIT
  all_paths = ["HLT_Mu7_IP4", "HLT_Mu8_IP6", "HLT_Mu8_IP5" , "HLT_Mu8_IP3", "HLT_Mu8p5_IP3p5", "HLT_Mu9_IP6", "HLT_Mu9_IP5", "HLT_Mu9_IP4", "HLT_Mu10p5_IP3p5", "HLT_Mu12_IP6"]
  #all_paths = ["HLT_Mu7_IP4", "HLT_Mu8_IP6"]
  ####

  if opt.whichana=='lumiTot' or opt.whichana=='lumiTotByPart' or opt.whichana=='lumiEach':
    paths_to_use = all_paths

  if opt.whichana=='lumiOverlap':
    #paths_to_use =  ["HLT_Mu7_IP4","HLT_Mu8_IP6"]   
    #paths_to_use =  ['HLT_Mu7_IP4']   
    paths_to_use = all_paths

  #### DO NOT EDIT
  dobyls = True          # recommended, do not change
  doCheckByTime = False if opt.whichana == 'overlaplumi' else True
  # if false, checks by run, fill and lumi_section,  check by time is faster and it gives same results as false if dobyls==True
  ####

  if opt.whichana=='lumiTotByPart' and opt.part==None: raise RuntimeError('with analysis lumiTotByPart enabled, you need to specify which part')

  if opt.whichana=='lumiTotByPart':
    fout = open('analysis_output_time_byls_lumiByPart_{part}.txt'.format(part=opt.part), 'w')
  else:
    fout = open('analysis_output_time_byls_{wa}.txt'.format(wa=opt.whichana), 'w')

  if dobyls:
    if opt.part != None:
      fname = 'briloutputs/output_byls_{part}'.format(part=opt.part)
    else:
      fname = 'briloutputs/output_byls'
    fname += '_{path}.csv'
  else:
    fname = 'briloutputs/output_{path}.csv'

  #### DO NOT EDIT
  periods = ['A', 'B', 'C', 'D']
  period_runs = {}
  period_runs['A'] = [315974, 316995]
  period_runs['B'] = [317087, 318877]
  period_runs['C'] = [319337, 320065]
  period_runs['D'] = [320674, 325172]
  ####

  ### create the dataframes
  dfs = {}
  print '==> Will read data from'
  for path in paths_to_use:
    df = pd.read_csv(fname.format(path=path), sep=',', comment='#') 
    print fname.format(path=path)
    dfs[path]=df
  # same calculation but doing the sum, so that it can be done for different periods


  ##########################
  # LumiEach (lumi for each trigger path)
  #########################
  if opt.whichana=='lumiEach':
    ## main analysis
    print '\n===> Lumi Each analysis'
    fout.write('path period lumi')

    lumi_each_period = {} 
    for path in paths_to_use:
      print '\n=====> path =', path 
      lumi_each_period[path] = {}
      lumi_each_period[path]['sum'] = 0
      for period in periods:
        lumi_each_period[path][period] = 0
      #
      this_df = dfs[path]
      for index, row in this_df.iterrows():
        #print row['run:fill']
        this_run = int( row['run:fill'].split(':')[0] )
        lumi_each_period[path]['sum'] += row['recorded(/fb)'] 
        for period in periods:
          if this_run >= period_runs[period][0] and this_run <= period_runs[period][1]:
            lumi_each_period[path][period] += row['recorded(/fb)']
      
      for period in periods:
        print '\n      period={}, lumi={:.3f}'.format(period,lumi_each_period[path][period]) 
      print '\n      period={}, lumi={:.3f}'.format('sum',lumi_each_period[path]['sum']) 
 
    for path in paths_to_use:
      for period in periods:
        fout.write('\n{:20s} {:10s} {:20.3f}'.format(path,period,lumi_each_period[path][period]))
      fout.write('\n{:20s} {:10s} {:20.3f}'.format(path,'sum',lumi_each_period[path]['sum']))

    ## cross-check calculation
    print '\n===> Cross-check analysis for all lumi analysis'
    lumis = {}
    for path in all_paths:
      with open(fname.format(path=path), 'r') as infile:
        for line in csv.reader(infile):
          # get the total
          if 'Sum recorded' in line[0]: 
            lumi=float(line[0].split(' : ')[1])
            #lumis.append(lumi)
            lumis[path] = {}
            lumis[path][path]=lumi
            to_print = '{:20s}  {:7.3f}'.format(path,lumi)
            print to_print
            ####fout.write(to_print + '\n')
  
  #########################
  # overlap lumis   
  #########################I
  if opt.whichana == 'lumiOverlap':

    lumi_overlap = {}
    print '\n===> Overlap Lumi Analysis' 
    for this_path in paths_to_use:
      print '\n =====> Check overlap of ', this_path
      this_df = dfs[this_path]
      lumi_overlap[this_path] = {}
  
      for check_path in paths_to_use:
        lumi_overlap[this_path][check_path] = 0.
        if check_path == this_path: continue
        print '        with:' , check_path
        check_df = dfs[check_path]
  
        # loop over rows of this_df and look for occurrences in the check_df
        past_times = []
        
        counter = 0
        for index, row in this_df.iterrows():
          counter +=1
          #if counter >= 5: break
          #this_fill = row['run:fill']
          this_time = row['time']
          if this_time in past_times: continue # do avoid double counting fills associated with different "parts"
          #print '  past_times', past_times
          #match_df = check_df.loc[check_df['run:fill']==this_fill & check_df['time']==this_time]
          #match_df = check_df.loc[check_df['run:fill']==this_fill]
          match_df = check_df.loc[check_df['time']==this_time]
          #print '  match_df' , match_df
          match_lumi = match_df['recorded(/fb)'].sum() if not match_df.empty else 0.
          #print '  match_lumi' , match_lumi
          #if not match_df.empty:
          #  print '  this_fill' , this_fill
          #  print '  this_time' , this_time
          #  print '  match_df' , match_df
  
          lumi_overlap[this_path][check_path] +=match_lumi # the sum goes over the fills associated to this_df
          #print '  tot_match_lumi' , lumi_overlap[this_path][check_path]
          past_times.append(this_time)
  
    #print ''
    print '          overlap lumi = ', lumi_overlap[this_path][check_path]
    #print '' 
  
    ## nice format printing
    header = []
    rows = []
    for p1 in paths_to_use:
      header.append(p1.split('HLT_')[1])
      row = []
      #row.append(p1.split('HLT_')[1])
      for p2 in paths_to_use:
        row.append('{:7.3f}'.format(lumi_overlap[p1][p2]))
      rows.append(row)
 
    print '===> Results' 
    print (' '.join(header))
    fout.write(' '.join(header))
    fout.write('\n')
    for row in rows:
      print( ' '.join(row))
      fout.write(' '.join(row))

  ########
  ## consistency check
  #
  #for path in paths_to_use:
  #  this_df = dfs[path]
  #  print '\n', path, 'tot lumi from sum : {:.3f}'.format(this_df['recorded(/fb)'].sum())



  #########################
  #### total lumi 
  #########################
  if (opt.whichana=='lumiTot' or opt.whichana=='lumiTotByPart'):
    # loop over all output csvs, sum lumi only if the lumiblock/time was not already counted from a previous trigger
    print '\n===> LumiTot analysis'
    if opt.whichana=='lumiTotByPart': print '\n     for {}'.format(opt.part)
    
    times_counted = {}
    run_fill_ls_counted = {}
    
    lumi_tot = 0.
    
    lumi_tot_period = {}
    lumi_tot_period['A'] = 0
    lumi_tot_period['B'] = 0
    lumi_tot_period['C'] = 0
    lumi_tot_period['D'] = 0
 
    for this_path in paths_to_use:
      print '\n====> path={}'.format(this_path)
      fout.write('\n' + this_path + '\n')
      this_df = dfs[this_path]
      times_counted[this_path] = []
      run_fill_ls_counted[this_path] = []
 
     
      for index, row in this_df.iterrows():    

        if doCheckByTime:
          this_time = row['time']
          skip_this_time = False
          # check if this_time is to be kept
          #skip_this_time = True
          for key in times_counted.keys():
            if this_path != key: 
              if this_time in times_counted[key]:
                skip_this_time = True
                break 
    
          if not skip_this_time:
            lumi_tot += row['recorded(/fb)']
            times_counted[this_path].append(this_time)
  
            # fill in also the lumi by period
            this_run = int( row['run:fill'].split(':')[0] )
            for period in periods:
              if this_run >= period_runs[period][0] and this_run <= period_runs[period][1]:
                lumi_tot_period[period] += row['recorded(/fb)']

        else: # check by lumi section
          this_run_fill = row['run:fill']
          this_ls = row['ls']
          this_run_fill_ls = this_run_fill + ' ' + this_ls
          #print(this_run_fill_ls)
          skip_this_run_fill_ls = False
          for key in run_fill_ls_counted.keys():
            if this_path != key: 
              #print(this_path, key)
              if this_run_fill_ls in run_fill_ls_counted[key]:
                skip_this_run_fill_ls = True
                #print(this_run_fill_ls, ' will be skipped ')
                break 
    
          if not skip_this_run_fill_ls:
            lumi_tot += row['recorded(/fb)']
            run_fill_ls_counted[this_path].append(this_run_fill_ls)
            
            # fill in also the lumi by period
            this_run = int( row['run:fill'].split(':')[0] )
            for period in periods:
              if this_run >= period_runs[period][0] and this_run <= period_runs[period][1]:
                lumi_tot_period[period] += row['recorded(/fb)']
     
      print '      lumi so far =' , lumi_tot
      print '      lumi by period so far:'
      for p,l in lumi_tot_period.items():
        print '         ',p,l

      fout.write('lumi so far =' + str(lumi_tot) + '\n')

    print '\nTotal Lumi for the B-parking dataset = {:.3f}'.format(lumi_tot)
    print '\nLumi by period'
    for p,l in lumi_tot_period.items():
      print '{}: {:.3f}'.format(p,l)

    fout.write('\nRESULTS:')
    fout.write('\nTotal Lumi for the B-parking dataset = {:.3f}\n'.format(lumi_tot))
    fout.write('\nLumi by period')
    for p,l in lumi_tot_period.items():
      fout.write('\n{}: {:.3f}'.format(p,l))
  fout.close()










