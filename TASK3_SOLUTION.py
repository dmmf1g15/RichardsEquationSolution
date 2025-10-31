#############
#Choose soil 
p=vg.GuelphLoamDrying
#Choose soil depth
ProfileDepth=3 #m

#Choose locations of soil in order to get weather. Ideally this location should have the soil-properties defined above. 
lat=55.917162
long=-3.187954

#Choose how long to do the simulation
days_in_past=50

# Initial conditions and make soil
dz=0.1
z=np.arange(dz/2.0,ProfileDepth,dz)
n=z.size
psi0=-z

# Grid in time
t = np.linspace(0,days_in_past,101) #solve for 30 days, saving the solution 101 times
dt=t[1]-t[0] #time step

#Make transient infiltration flux
t_interp=np.arange(0,t[-1]+2,1) #need an extra day for the rainfall
dates,netflux=get_historic_netflux(days_in_past+2,-1,len(t_interp))
qTfun=interp1d(t_interp,-1*netflux) #-1 for flux condition since we are in upside down land

#Make plot showing rainfall profile and make sure interpolation agrees
plt.figure(figsize=(10, 5))
plt.plot(t,-qTfun(t))
plt.scatter(dates,netflux,color='r')
plt.ylabel('ET [m/day]')
plt.xticks(rotation=45, ha='right',fontsize=6)
plt.xlabel('day')


##RUN THE MODEL
psi=Solve_Richards(psi0,t,dz,n,p,qTfun)


###############
###Plots

#Post process model output to get useful information
# Get water content from pressure
theta=vg.thetaFun(psi,p)

# Get total profile storage
S=theta.sum(axis=1)*dz

# Get change in storage [dVol]
dS=np.zeros(S.size)
dS[1:]=np.diff(S)/(t[1]-t[0])

#rainfall
qI=qTfun(t)

# Get discharge flux
KBot=vg.KFun(psi[:,0],p)
qD=-KBot


plt.figure(figsize=(10, 10))
#Make mass balance plots
plt.subplot(311)
plt.plot(t,S)
plt.xlabel("day")
plt.ylabel("Volume of water in soil [m]")
plt.subplot(312)
plt.plot(t,dS,label='Change in storage')
plt.plot(t,-qI,label='Infiltration')
#plt.plot(t,-qD,label='Discharge')
plt.ylabel('m/day')
plt.xlabel('day')
plt.legend()
