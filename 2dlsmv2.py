
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import numpy as np 
import time
option = st.sidebar.selectbox("Which Dashboard?", ('Home','Inhomogeneity Example', 'Input Morphology'), 0)
st.header(option)
if option == 'Home':
	st.subheader('2D LSM Informaiton')
	st.write('The Inhomogeneity Example page will create a model that will capture the elastic deformation of a particle in a matrix.')
	st.write('If you already have have a data set ready then you can use the Input Morphology to upload your Excel, CSV, or .txt file. To get the program to run correctly make sure it is formatted the same as the data file from the model page. It should have no column headers. The first column should contain x cordinates, second should contain y cordinates and the third should contain elastic. The radius is set to 10 on this model.')

if option == 'Inhomogeneity Example':

	# A regular 2D Lattice Spring Model. 100x100 and using conjugate gradient solver.
	SIZE = 100

	# The value of residual forces when the system has reached a solution
	EPSILON = 1e-5

	# The force applied to the boundary nodes (negative to left and positive to right)
	FORCE = 0.1
	
	# The elastic moduli - from 1 to 10. Assuming C's = 0 for now
	userELASTIC = st.slider("Enter a value for Elastic ",  min_value=0.0, max_value=10.0, value=1.0, step=.5)
	ELASTIC = [[1.0 for col in range(SIZE) ] for row in range(SIZE)]
	num = st.slider("Enter a value for radius ", min_value=10, max_value=20, value=10, step=5)
	if num < 10:
		rad = 100.
	elif num ==15:
		rad =225.
	else:
		rad = 300.
	#if (X - SIZE/2)*(X-SIZE/2) + (Y-SIZE/2)*(Y-SIZE/2) < R*R:
	#elastic = elastic+1
	#The solution requires the displacements and other variables
	DX = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	DY = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	XR = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	YR = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	XAP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	YAP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	XP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	YP = [[0. for col in range(SIZE) ] for row in range(SIZE)]

	IX = [0 for col in range(9) ]
	IY = [0 for col in range(9) ]
	
	STRAIN = [[0. for col in range(SIZE) ] for row in range(SIZE)]
	
	# This is where the Elastic constants have to be loaded in. Say, from 1 to 10 in terms of the values. 
	X = 0
	while(X < SIZE):
		Y = 0
		while(Y < SIZE):
			ELASTIC[X][Y] = 1.
			R = (X-0.5*SIZE)*(X-0.5*SIZE)+(Y-0.5*SIZE)*(Y-0.5*SIZE)
			if(R < rad): 
				ELASTIC[X][Y] = userELASTIC
			Y = Y + 1
		X = X + 1

	# Need to find neighbors
	IX[1] = 1
	IY[1] = 0
	IX[2] = -1
	IY[2] = 0
	IX[3] = 0
	IY[3] = 1
	IX[4] = 0
	IY[4] = -1

	IX[5] = 1
	IY[5] = 1
	IX[6] = -1
	IY[6] = -1
	IX[7] = 1
	IY[7] = -1
	IX[8] = -1
	IY[8] = 1
	# Calculate Initial R
	if st.button("Run"):
		with st.spinner('Calculating Results...'):
			RR = 0.
			X = 0
			while(X < SIZE):
				Y = 0
				while(Y < SIZE):
					
					AX = 0.
					AY = 0.
					FX = 0.
					FY = 0.
					
					i = 1
					while(i < 9):
						X2 = X + IX[i]
						Y2 = Y + IY[i]
						
						if(X2 < 0 or X2 > SIZE-1): 
							i = i + 1
							continue
						
						if(Y2 < 0 or Y2 > SIZE-1): 
							i = i + 1
							continue
					
						K = 0.5*(ELASTIC[X][Y] + ELASTIC[X2][Y2])
						
						if(K < 1.): K = 1.
						if(K > 10.): K = 10.
						
						if(i > 4):
							if(X == 0 and Y == 0): 
								K = K * 2.
							if(X == SIZE-1 and Y == 0): 
								K = K * 2.
							if(X == 0 and Y == SIZE-1): 
								K = K * 2.
							if(X == SIZE-1 and Y == SIZE-1): 
								K = K * 2.
							
						dx = DX[X2][Y2] - DX[X][Y]
						dy = DY[X2][Y2] - DY[X][Y]

						if(i == 1 or i == 2): 
							AX -= K * dx

						if(i == 3 or i == 4): 
							AY -= K * dy

						if(i == 5 or i == 6): 
							AX -= 0.25 * (K * dx + K * dy)
							AY -= 0.25 * (K * dx + K * dy)

						if(i == 7 or i == 8): 
							AX -= 0.25 * (K * dx - K * dy)
							AY -= 0.25 * (-K * dx + K * dy)
										
						i = i + 1
						
					if(X == 0): 
						FX = -FORCE
					
					if(X == SIZE-1): 
						FX = FORCE
					
					XR[X][Y] = XP[X][Y] = FX - AX
					YR[X][Y] = YP[X][Y] = FY - AY
					
					RR += XR[X][Y] * XR[X][Y]
					RR += YR[X][Y] * YR[X][Y]
					
					Y = Y + 1
				X = X + 1

			print('Residual = ', format(RR,'.3E'), end="\r")

			MIN = RR
			# Iteration until desired precision is obtained
			with st.empty():
				while(RR > EPSILON):
					
					PAP = 0.
					
					X = 0
					while(X < SIZE):
						Y = 0
						while(Y < SIZE):
						
							XAP[X][Y] = 0.
							YAP[X][Y] = 0.
						
							i = 1
							while(i < 9):
								X2 = X + IX[i]
								Y2 = Y + IY[i]
							
								if(X2 < 0 or X2 > SIZE-1): 
									i = i + 1
									continue
							
								if(Y2 < 0 or Y2 > SIZE-1): 
									i = i + 1
									continue
						
								K = 0.5*(ELASTIC[X][Y] + ELASTIC[X2][Y2])
							
								if(K < 1.): K = 1.
								if(K > 10.): K = 10.
							
								if(i > 4):
									if(X == 0 and Y == 0): 
										K = K * 2.
									if(X == SIZE-1 and Y == 0): 
										K = K * 2.
									if(X == 0 and Y == SIZE-1): 
										K = K * 2.
									if(X == SIZE-1 and Y == SIZE-1): 
										K = K * 2.
								
								dx = XP[X2][Y2] - XP[X][Y]
								dy = YP[X2][Y2] - YP[X][Y]

								if(i == 1 or i == 2): 
									XAP[X][Y] -= K * dx

								if(i == 3 or i == 4): 
									YAP[X][Y] -= K * dy

								if(i == 5 or i == 6): 
									XAP[X][Y] -= 0.25 * (K * dx + K * dy)
									YAP[X][Y] -= 0.25 * (K * dx + K * dy)

								if(i == 7 or i == 8): 
									XAP[X][Y] -= 0.25 * (K * dx - K * dy)
									YAP[X][Y] -= 0.25 * (-K * dx + K * dy)
											
								i = i + 1
							
							PAP += XP[X][Y] * XAP[X][Y]
							PAP += YP[X][Y] * YAP[X][Y]
						
							Y = Y + 1
						X = X + 1
						
					ALPHA = -RR/PAP
					
					RR1 = 0.0
					
					X = 0
					while(X < SIZE):
						Y = 0
						while(Y < SIZE):
						
							DX[X][Y] -= ALPHA * XP[X][Y]
							DY[X][Y] -= ALPHA * YP[X][Y]
					
							XR[X][Y] += ALPHA * XAP[X][Y]
							YR[X][Y] += ALPHA * YAP[X][Y]
							
							RR1 += XR[X][Y] * XR[X][Y]
							RR1 += YR[X][Y] * YR[X][Y]
					
							Y = Y + 1
						X = X + 1

					st.write('Residual = ', format(RR1,'.3E'), 'Minimum = ', format(MIN,'.3E'), '          ', end="\r")

					if(RR1 < MIN): MIN = RR1
					
					BETA = RR1/RR
					
					X = 0
					while(X < SIZE):
						Y = 0
						while(Y < SIZE):
						
							XP[X][Y] = XR[X][Y] + BETA * XP[X][Y]
							YP[X][Y] = YR[X][Y] + BETA * YP[X][Y]
					
							Y = Y + 1
						X = X + 1
					
					RR = RR1
					
					if(RR > 1000.*MIN):
						print('Looks like system is not converging')
					

				# Once it's solved we have to decide what to output!

				# Young's Modulus and Poissons ratio
				XSYST = 0.
				YSYST = 0.
					
				X = 0
				while(X < SIZE):
					Y = 0
					while(Y < SIZE):
						if(X == 0): XSYST -= DX[X][Y]
						if(X == SIZE-1): XSYST += DX[X][Y]
						if(Y == 0): YSYST -= DY[X][Y]
						if(Y == SIZE-1): YSYST += DY[X][Y]

						Y = Y + 1
					X = X + 1

				XSYST /= SIZE
				YSYST /= SIZE

				YOUNG = (SIZE - 1.) * (FORCE / XSYST)
				POISSON = YSYST/XSYST
		st.success('Done!')
		#report section 
		st.write('Residual = ', RR)
		st.write('Youngs modulus = ', YOUNG)
		st.write('Poissons ratio = ', POISSON)
		print('\n')

		# Calculate strain in tensile direction only 
		X = 0
		while(X < SIZE):
			Y = 0
			while(Y < SIZE):
				if(X < 1): STRAIN[X][Y] = DX[X+1][Y] - DX[X][Y]
				elif(X > SIZE -2): STRAIN[X][Y] = DX[X][Y] - DX[X-1][Y]
				else: STRAIN[X][Y] = 0.5*(DX[X+1][Y] - DX[X-1][Y])
					
				Y = Y + 1
			X = X + 1

		OUTPUT = open("lsm.txt", "w")
	
		X = 0
		while(X < SIZE):
			Y = 0
			while(Y < SIZE):
				string = str(X) + ' ' + str(Y) + ' ' + str(STRAIN[X][Y]) + "\n"
				OUTPUT.write(string)
				Y = Y + 1
			OUTPUT.write(string)
			X = X + 1
		
		OUTPUT = open("elasticlsm.txt", "w")
	
		X = 0
		while(X < SIZE):
			Y = 0
			while(Y < SIZE):
				string = str(X) + ' ' + str(Y) + ' ' + str(ELASTIC[X][Y]) + "\n"
				OUTPUT.write(string)
				Y = Y + 1
			OUTPUT.write(string)
			X = X + 1
		
		pts = np.loadtxt('lsm.txt',dtype=float, delimiter=' ')
		x, y, z = pts.T
		

		fig2 = go.Figure(data = go.Contour(z=z, x=x, y=y))
		st.write(fig2)

		with open('lsm.txt') as f:
			st.download_button('Download strain data file', f,'result.csv', 'text/csv')
		with open('elasticlsm.txt') as f:
			st.download_button('Download elastic data file', f,'elasticresult.csv', 'text/csv')


if option == 'Input Morphology':
	uploaded_file =  st.sidebar.file_uploader(label="upload your excel file here.", type =['xlsx','csv','txt'])
	if uploaded_file is not None:
		try:
			df = pd.read_csv(uploaded_file, sep=" ", header=None)
			df.to_csv('results2d.txt',index=False, sep=" ")
			pts = np.loadtxt('results2d.txt',dtype=float, delimiter=' ')
			x1, y1, ela = pts.T
			if df is not None:
				st.write(df)
				# A regular 2D Lattice Spring Model. 100x100 and using conjugate gradient solver.
				SIZE = 100

				# The value of residual forces when the system has reached a solution
				EPSILON = 1e-5

				# The force applied to the boundary nodes (negative to left and positive to right)
				FORCE = 0.1
				
				# The elastic moduli - from 1 to 10. Assuming C's = 0 for now
				# The elastic moduli - from 1 to 10. Assuming C's = 0 for now
				l = list(ela)
				ELASTIC = [[1. for col in range(SIZE) ] for row in range(SIZE)]
				i = 0
				for n in l:
					X = i % SIZE
					Y = i % SIZE
					ELASTIC[X][Y] = n                       
					i=i+1
				 
				#The solution requires the displacements and other variables
				DX = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				DY = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				XR = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				YR = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				XAP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				YAP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				XP = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				YP = [[0. for col in range(SIZE) ] for row in range(SIZE)]

				IX = [0 for col in range(9) ]
				IY = [0 for col in range(9) ]
				STRAIN = [[0. for col in range(SIZE) ] for row in range(SIZE)]
				
				# Need to find neighbors
				IX[1] = 1
				IY[1] = 0
				IX[2] = -1
				IY[2] = 0
				IX[3] = 0
				IY[3] = 1
				IX[4] = 0
				IY[4] = -1

				IX[5] = 1
				IY[5] = 1
				IX[6] = -1
				IY[6] = -1
				IX[7] = 1
				IY[7] = -1
				IX[8] = -1
				IY[8] = 1
				#st.write(ELASTIC)
				# Calculate Initial R
				if st.button("Run"):
					with st.spinner('Calculating Results...'):
						RR = 0.
						X = 0
						while(X < SIZE):
							Y = 0
							while(Y < SIZE):
								
								AX = 0.
								AY = 0.
								FX = 0.
								FY = 0.
								
								i = 1
								while(i < 9):
									X2 = X + IX[i]
									Y2 = Y + IY[i]
									
									if(X2 < 0 or X2 > SIZE-1): 
										i = i + 1
										continue
									
									if(Y2 < 0 or Y2 > SIZE-1): 
										i = i + 1
										continue
								
									K = 0.5*(ELASTIC[X][Y] + ELASTIC[X2][Y2])
									
									if(K < 1.): K = 1.
									if(K > 10.): K = 10.
									
									if(i > 4):
										if(X == 0 and Y == 0): 
											K = K * 2.
										if(X == SIZE-1 and Y == 0): 
											K = K * 2.
										if(X == 0 and Y == SIZE-1): 
											K = K * 2.
										if(X == SIZE-1 and Y == SIZE-1): 
											K = K * 2.
										
									dx = DX[X2][Y2] - DX[X][Y]
									dy = DY[X2][Y2] - DY[X][Y]

									if(i == 1 or i == 2): 
										AX -= K * dx

									if(i == 3 or i == 4): 
										AY -= K * dy

									if(i == 5 or i == 6): 
										AX -= 0.25 * (K * dx + K * dy)
										AY -= 0.25 * (K * dx + K * dy)

									if(i == 7 or i == 8): 
										AX -= 0.25 * (K * dx - K * dy)
										AY -= 0.25 * (-K * dx + K * dy)
													
									i = i + 1
									
								if(X == 0): 
									FX = -FORCE
								
								if(X == SIZE-1): 
									FX = FORCE
								
								XR[X][Y] = XP[X][Y] = FX - AX
								YR[X][Y] = YP[X][Y] = FY - AY
								
								RR += XR[X][Y] * XR[X][Y]
								RR += YR[X][Y] * YR[X][Y]
								
								Y = Y + 1
							X = X + 1

						print('Residual = ', format(RR,'.3E'), end="\r")

						MIN = RR
						# Iteration until desired precision is obtained
						with st.empty():
							while(RR > EPSILON):
								
								PAP = 0.
								
								X = 0
								while(X < SIZE):
									Y = 0
									while(Y < SIZE):
									
										XAP[X][Y] = 0.
										YAP[X][Y] = 0.
									
										i = 1
										while(i < 9):
											X2 = X + IX[i]
											Y2 = Y + IY[i]
										
											if(X2 < 0 or X2 > SIZE-1): 
												i = i + 1
												continue
										
											if(Y2 < 0 or Y2 > SIZE-1): 
												i = i + 1
												continue
									
											K = 0.5*(ELASTIC[X][Y] + ELASTIC[X2][Y2])
										
											if(K < 1.): K = 1.
											if(K > 10.): K = 10.
										
											if(i > 4):
												if(X == 0 and Y == 0): 
													K = K * 2.
												if(X == SIZE-1 and Y == 0): 
													K = K * 2.
												if(X == 0 and Y == SIZE-1): 
													K = K * 2.
												if(X == SIZE-1 and Y == SIZE-1): 
													K = K * 2.
											
											dx = XP[X2][Y2] - XP[X][Y]
											dy = YP[X2][Y2] - YP[X][Y]

											if(i == 1 or i == 2): 
												XAP[X][Y] -= K * dx

											if(i == 3 or i == 4): 
												YAP[X][Y] -= K * dy

											if(i == 5 or i == 6): 
												XAP[X][Y] -= 0.25 * (K * dx + K * dy)
												YAP[X][Y] -= 0.25 * (K * dx + K * dy)

											if(i == 7 or i == 8): 
												XAP[X][Y] -= 0.25 * (K * dx - K * dy)
												YAP[X][Y] -= 0.25 * (-K * dx + K * dy)
														
											i = i + 1
										
										PAP += XP[X][Y] * XAP[X][Y]
										PAP += YP[X][Y] * YAP[X][Y]
									
										Y = Y + 1
									X = X + 1
									
								ALPHA = -RR/PAP
								
								RR1 = 0.0
								
								X = 0
								while(X < SIZE):
									Y = 0
									while(Y < SIZE):
									
										DX[X][Y] -= ALPHA * XP[X][Y]
										DY[X][Y] -= ALPHA * YP[X][Y]
								
										XR[X][Y] += ALPHA * XAP[X][Y]
										YR[X][Y] += ALPHA * YAP[X][Y]
										
										RR1 += XR[X][Y] * XR[X][Y]
										RR1 += YR[X][Y] * YR[X][Y]
								
										Y = Y + 1
									X = X + 1

								st.write('Residual = ', format(RR1,'.3E'), 'Minimum = ', format(MIN,'.3E'), '          ', end="\r")

								if(RR1 < MIN): MIN = RR1
								
								BETA = RR1/RR
								
								X = 0
								while(X < SIZE):
									Y = 0
									while(Y < SIZE):
									
										XP[X][Y] = XR[X][Y] + BETA * XP[X][Y]
										YP[X][Y] = YR[X][Y] + BETA * YP[X][Y]
								
										Y = Y + 1
									X = X + 1
								
								RR = RR1
								
								if(RR > 1000.*MIN):
									print('Looks like system is not converging')
								

							# Once it's solved we have to decide what to output!

							# Young's Modulus and Poissons ratio
							XSYST = 0.
							YSYST = 0.
								
							X = 0
							while(X < SIZE):
								Y = 0
								while(Y < SIZE):
									if(X == 0): XSYST -= DX[X][Y]
									if(X == SIZE-1): XSYST += DX[X][Y]
									if(Y == 0): YSYST -= DY[X][Y]
									if(Y == SIZE-1): YSYST += DY[X][Y]

									Y = Y + 1
								X = X + 1

							XSYST /= SIZE
							YSYST /= SIZE

							YOUNG = (SIZE - 1.) * (FORCE / XSYST)
							POISSON = YSYST/XSYST
					st.success('Done!')
					#report section 
					st.write('Residual = ', RR)
					st.write('Youngs modulus = ', YOUNG)
					st.write('Poissons ratio = ', POISSON)
					print('\n')

					# Calculate strain in tensile direction only 
					X = 0
					while(X < SIZE):
						Y = 0
						while(Y < SIZE):
							if(X < 1): STRAIN[X][Y] = DX[X+1][Y] - DX[X][Y]
							elif(X > SIZE -2): STRAIN[X][Y] = DX[X][Y] - DX[X-1][Y]
							else: STRAIN[X][Y] = 0.5*(DX[X+1][Y] - DX[X-1][Y])
								
							Y = Y + 1
						X = X + 1

					OUTPUT = open("lsm2.txt", "w")
				
					X = 0
					while(X < SIZE):
						Y = 0
						while(Y < SIZE):
							string = str(X) + ' ' + str(Y) + ' ' + str(STRAIN[X][Y]) + "\n"
							OUTPUT.write(string)
							Y = Y + 1
						OUTPUT.write(string)
						X = X + 1
					
					pts = np.loadtxt('lsm2.txt',dtype=float, delimiter=' ')
					x, y, z = pts.T
					

					fig3 = go.Figure(data = go.Contour(z=z, x=x, y=y))
					st.write(fig3)
					with open('lsm2.txt') as f:
						st.download_button('Download dat file', f,'result.csv', 'text/csv')
		except Exception as e:
			st.write("file not in correct file please check example page for file format.")
			#df = pd.read_csv(uploaded_file)
			
	















