#!/usr/bin/env python3
# import time
import numpy as np
import matplotlib.pyplot as plt

#from Landlab
from landlab import RasterModelGrid, imshow_grid, imshowhs_grid
from landlab.io import read_esri_ascii

#Hillslope geomorphology
from landlab.components import ExponentialWeatherer
from landlab.components import DepthDependentTaylorDiffuser
from landlab.components import DepthDependentDiffuser

#Fluvial Geomorphology and Flow routing
from landlab.components import FlowDirectorMFD #trying the FlowDirectorMFD
from landlab.components import FlowAccumulator, Space, FastscapeEroder, PriorityFloodFlowRouter
from landlab.components.space import SpaceLargeScaleEroder

from ss_fault_function import ss_fault

#READING STEADY STATE TOPO
(mg,z)=read_esri_ascii('/Users/taranguiz/Research/CSDMS_summer_2022/output_new_topo_ddd_5/finaltopo_topographic__elevation.asc', name="topographic__elevation")
(mg1,soil_0)=read_esri_ascii('/Users/taranguiz/Research/CSDMS_summer_2022/output_new_topo_ddd_5/finaltopo_soil__depth.asc', name='soil__depth')
(mg2,bed_0)=read_esri_ascii('/Users/taranguiz/Research/CSDMS_summer_2022/output_new_topo_ddd_5/finaltopo_bedrock__elevation.asc', name='bedrock__elevation')
# (mg,soil_prod)=read_esri_ascii('finaltopo_soil_production__rate.asc', name='soil_production__rate')
#print(mg.at_node.keys())

mg.add_field("soil__depth", soil_0, at='node')
mg.add_field("bedrock__elevation", bed_0, at='node')


mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=False, left_is_closed=True,
                                       right_is_closed=True, top_is_closed=True)
shrink = 0.5
imshow_grid(mg,z, cmap='viridis', shrink=shrink)
plt.title('Initial Topography')
plt.show()

#according to info from the file
ymax=1000
xmax=3000
dxy=10 #spacing

nrows = int(ymax/dxy)
ncols = int(xmax/dxy)


# mg.add_zeros("node", "soil__depth") #add field to the grid
# mg.at_node["soil__depth"]=mg.at_node["soil__depth"]+2 #2
# mg.at_node["bedrock__elevation"]=mg.at_node["topographic__elevation"] - mg.at_node["soil__depth"]

soil=mg.at_node['soil__depth']
bed= mg.at_node['bedrock__elevation']

print (mg.at_node['topographic__elevation'][mg.core_nodes])
print (mg.at_node['bedrock__elevation'][mg.core_nodes])
print (mg.at_node['soil__depth'][mg.core_nodes])

# Geomorphic parameters
# uplift
uplift_rate= 3*1e-5

#Hillsope Geomorphology for DDTD component
# H=5 # original soil depth
Sc=0.7
Hstar= 0.1 # characteristic transport depth, m
V0= 0.1 #transport velocity coefficient
D= Hstar*V0#V0 *Hstar  #effective(maximum) diffusivity

#Fluvial Erosion for SPACE Large Scale Eroder
K_sed=5*1e-5 #sediment erodibility
K_br= 1*1e-5 #bedrock erodibility
F_f=0.5 #fraction of fine sediment
phi= 0.5 #sediment porosity
H_star=Hstar #sediment entrainment lenght scale
Vs= 1 #velocity of sediment
m_sp= 0.3 #exponent ondrainage area stream power
n_sp= 1 #exponent on channel slope in the stream power framework
sp_crit_sed=0 #sediment erosion threshold
sp_crit_br=0 #bedrock erosion threshold

#instantiate components
expweath=ExponentialWeatherer(mg, soil_production__maximum_rate=1*1e-5, soil_production__decay_depth=Hstar)

# Hillslope with Taylor Diffuser
ddtd=DepthDependentTaylorDiffuser(mg,slope_crit=Sc,
                                  soil_transport_velocity=V0,
                                  soil_transport_decay_depth=Hstar,
                                  nterms=2,
                                  dynamic_dt=True,
                                  if_unstable='warn')
#Flow Router
fr=PriorityFloodFlowRouter(mg, flow_metric='D8', suppress_out=True, runoff_rate=0.5)
#SPACE Large Scale
space= SpaceLargeScaleEroder(mg,
                             K_sed=K_sed,
                             K_br=K_br,
                            F_f=F_f,
                            phi=phi,
                            H_star=Hstar,
                            v_s=Vs,
                            m_sp=m_sp,
                            n_sp=n_sp,
                            sp_crit_sed=0,
                             sp_crit_br=0)
#TECTONICS AND TIME PARAMETERS
total_slip=500
total_model_time= 500000
dt=100
iterations= np.arange(0,total_model_time,dt)
print(iterations)
desired_slip_per_event=(total_slip/total_model_time)*dt
shrink = 0.5
fault_loc_y=int(mg.number_of_node_rows / 3.)
fault_nodes = np.where(mg.node_y==(fault_loc_y*10))[0]
print(fault_nodes)
figsize = [16,4] # size of grid plots
fig, ax = plt.subplots(figsize=figsize)
x = mg.node_x[fault_nodes]
# surface_level = soil
# ax.plot(x, soil[fault_nodes], 'orange', linewidth=2, markersize=12, label='soil')
# ax.plot(x, bed[fault_nodes], linewidth=2, markersize=12, label='bedrock')
ax.plot(x, z[fault_nodes],'red', linewidth=2, markersize=12, label='topo')


plt.title('Original cross-Profile topography at fault location')
#plt.text(480, 9, 'air')
#plt.text(480, 7, 'soil')
#plt.text(470, 2, 'bedrock')

# plt.xlim(0,1000)
# plt.ylim(0, 10)
ax.set_xlabel('X (m)')
ax.set_ylabel('Depth (m)')
ax.legend(loc='lower right')
plt.show()

#fluvial array (look up table)
fluvial_freq=80000 #how often the humid period occurs
fluvial_len=5000 #how long the humid period last
fluvial_0=np.arange(fluvial_freq,total_model_time, fluvial_freq)
fluvial_n=fluvial_0+fluvial_len
fluvial_times=np.vstack((fluvial_0,fluvial_n)).T
print(fluvial_times)

#hillslope array
hillslope_freq=40000
hillslope_len=5000
hillslope_0=np.arange(hillslope_freq,total_model_time, hillslope_freq)
hillslope_n=hillslope_0+hillslope_len
hillslope_times=np.vstack((hillslope_0, hillslope_n)).T
print(hillslope_times)

#things for the loop
time=0 #time counter
f=0 #index counter for fluvial
h=0 #index counter for hillslope
accumulate=0

while time < total_model_time:

      z[mg.core_nodes]+= (uplift_rate*dt) #do uplift all the time
      bed[mg.core_nodes]+= (uplift_rate*dt)

      accumulate += desired_slip_per_event
      print('is accumulating')

      if accumulate >= mg.dx:
          ss_fault(grid=mg, fault_loc_y=fault_loc_y, total_slip=total_slip,
                   total_time=total_model_time, method='roll', accumulate=accumulate)
          accumulate = accumulate % mg.dx
          # expweath.maximum_weathering_rate=1*1e-6
          # expweath.calc_soil_prod_rate()
          # ddtd.run_one_step(dt=1250)
          print('one slip')

      if time >= fluvial_times[f,0] and time <= fluvial_times[f,1]:
          print('is: ' + str(time) +' so is hola fluvial')
          fr.run_one_step()
          space.run_one_step(dt)
          if time == fluvial_times[f,1] and f < (len(fluvial_times)-1):
              f=f+1
          if f==(len(fluvial_times)-1):
              pass
      if time >= hillslope_times[h,0] and time <= hillslope_times[h,1]:
          print ('is: ' + str(time) +' so is hola colluvial')
          expweath.maximum_weathering_rate=1*1e-5
          expweath.calc_soil_prod_rate()
          ddtd.run_one_step(dt)
          if time == hillslope_times[h,1] and h < (len(hillslope_times)-1):
              h=h+1
          if h == (len(hillslope_times)-1):
              pass

      # if time%5000 == 0:
      #     imshow_grid(mg, z, cmap='viridis', shrink=shrink)
      #     plt.title('Topography after ' + str(time+dt) + ' years')
      #     plt.show()

      time = time + dt
      print(time)

print(time)
imshow_grid(mg,z, cmap='viridis', shrink=shrink)
plt.title('Topography after ' + str(total_model_time) + ' years')
plt.show()

figsize = [16,4] # size of grid plots
fig, ax = plt.subplots(figsize=figsize)

x = mg.node_x[fault_nodes]
# soil_level = bed + soil

ax.plot(x, mg.at_node['soil__depth'][fault_nodes], 'orange', linewidth=2, markersize=12, label='soil')
ax.plot(x, mg.at_node['bedrock__elevation'][fault_nodes], linewidth=2, markersize=12, label='bedrock')
ax.plot(x, mg.at_node['topographic__elevation'][fault_nodes],'red', linewidth=2, markersize=12, label='topo')


plt.title('Final cross-Profile topography at fault location humid inducing hillslope activity with 1 mm/yr fault')
#plt.text(480, 9, 'air')
#plt.text(480, 7, 'soil')
#plt.text(470, 2, 'bedrock')

# plt.xlim(0,1000)
# plt.ylim(0, 10)
ax.set_xlabel('X (m)')
ax.set_ylabel('Depth (m)')
ax.legend(loc='upper right')
plt.show()









