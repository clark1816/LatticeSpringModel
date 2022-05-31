import plotly_express as px
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import numpy as np 
option = st.sidebar.selectbox("Which Dashboard?", ('home page', '2D Model Page', '2D Graph Page'), 2)
st.header(option)

if option == '2D Model Page':

    # A regular 2D Lattice Spring Model. 100x100 and using conjugate gradient solver.
    SIZE = 100

    # The value of residual forces when the system has reached a solution
    EPSILON = 1e-5

    # The force applied to the boundary nodes (negative to left and positive to right)
    FORCE = 0.1

    # The elastic moduli - from 1 to 10. Assuming C's = 0 for now
    ELASTIC = [[0. for col in range(SIZE) ] for row in range(SIZE)]

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
            #if(R < 100.): ELASTIC[X][Y] = 10.
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
            
            print('Residual = ', format(RR1,'.3E'), 'Minimum = ', format(MIN,'.3E'), '          ', end="\r")
            
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

    OUTPUT = open("lsm.dat", "w")
    X = 0
    #OUTPUT.write('x' + ',' + 'y' + ',' + 'z' + "\n")
    while(X < SIZE):
        Y = 0
        while(Y < SIZE):
            string = str(X) + ' ' + str(Y) + ' ' + str(STRAIN[X][Y]) + "\n"
            OUTPUT.write(string);
            Y = Y + 1
        OUTPUT.write(string);
        X = X + 1
    with open('lsm.dat') as f:
        st.download_button('Download dat file', f,'lsm2d.txt', 'text/csv')


if option == '2D Graph Page':
    st.write('upload excel file')
    uploaded_file =  st.sidebar.file_uploader(label="upload your excel file here.", type =['xlsx','csv','txt'])
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
        except Exception as e:
            print(e)
            df = pd.read_txt(uploaded_file)
    try:
        #plots mesh3d
        pts = np.loadtxt(np.DataSource().open('https://raw.githubusercontent.com/clark1816/LatticeSpringModel/main/2dlsm.txt'))
        x, y, z = pts.T

        fig = go.Figure(data=[go.Mesh3d(x=x, y=y, z=z,
                   alphahull=5,
                   opacity=0.4,
                   color='cyan')])
        fig.show()
        fig2 = go.Figure(data = 
            go.Contour(z=z, x=x, y=y))
        fig2.show()
    except Exception as e: 
        print(e)
        st.write('Please upload file to the application ')

# #Use matplotlib????
#  PX, PY = np.meshgrid(np.linspace(0, SIZE-1, SIZE), np.linspace(0, SIZE-1, SIZE))
#  Z = np.asarray(STRAIN)
#  levels = np.linspace(Z.min(), Z.max(), 7)

# # # plot
#  fig, ax = plt.subplots()
#  ax.contourf(PX, PY, Z, levels=levels)
#  plt.show()














