import numpy as np
import matplotlib.pyplot as plt
import math
from io import BytesIO
import streamlit as st
import pandas as pd
from matplotlib.backends.backend_agg import RendererAgg
_lock = RendererAgg.lock

of_dict = {
    'A'  : {'Name' : 'Organo Facies A - Oil' , 'A' : float(2.13e13) , 'E' : float(206.4 * 1e3 * 0.000239006) , 's' : float(8.2 * 1000 * 0.000239006)},
    'B'  : {'Name' : 'Organo Facies B - Oil' , 'A' : float(8.14e13) , 'E' : float(215.2 * 1e3 * 0.000239006) , 's' : float(8.3 * 1000 * 0.000239006)},
    'C'  : {'Name' : 'Organo Facies C - Oil' , 'A' : float(2.44e14) , 'E' : float(221.4 * 1e3 * 0.000239006) , 's' : float(3.9 * 1000 * 0.000239006)},
    'DE' : {'Name' : 'Organo Facies DE - Oil', 'A' : float(4.97e14) , 'E' : float(228.2 * 1e3 * 0.000239006) , 's' : float(7.9 * 1000 * 0.000239006)},
    'F'  : {'Name' : 'Organo Facies F - Oil' , 'A' : float(1.23e17) , 'E' : float(259.1 * 1e3 * 0.000239006) , 's' : float(6.6 * 1000 * 0.000239006)}
    }

def Arrhenius(A, Ea, T, time): # time is in seconds, T is in Kelvin
    return math.exp(-time * A * math.exp(-Ea / (0.001987 * (T + 273))))


def computeEzRo(vitrinitetble ,T, dt):
    for i in range(vitrinitetble.shape[0]):
        vitrinitetble[i][1] *= Arrhenius(1.0e13 , vitrinitetble[i][0] , T , dt)
    Tr = (85.0 - np.sum(vitrinitetble , axis = 0)[1]) / 100.0
    return math.exp(-1.6 + 3.7 * Tr) / 100 , vitrinitetble


def primarycracking(kerogen, alpha, user_data, top_fig_select, scalar_start, scalar_end):

    with _lock:
   
        labels = ["Temperature", "HI", "TR", "EasyRo", "TMax", "Thermal Stress", "Age"]
        np.set_printoptions(precision=2, suppress = True)
        temptable = np.linspace(30.,180.,16)
        vitrinitetble = np.array([[34.0 , 3.0],[36.0 , 3.0],[38.0 , 4.0],[40.0 , 4.0],[42.0 , 5.0],
                                [44.0 , 5.0],[46.0 , 6.0],[48.0 , 4.0],[50.0 , 4.0],[52.0 , 7.0],
                                [54.0 , 6.0],[56.0 , 6.0],[58.0 , 6.0],[60.0 , 5.0],[62.0 , 5.0],
                                [64.0 , 4.0],[66.0 , 3.0],[68.0 , 2.0],[70.0 , 2.0],[72.0 , 1.0]])
        T = 0
        T2 = 0
        age = 0
        sInMa = 60 * 60 * 24 *365 * 1e6
        result = []
        result_dict = {}
        maxT = 220
        dt = 1 / alpha #dt in Ma

        dts = sInMa * dt
        TMax = 0

        title = str(kerogen.name) + ' \n\n Kinetic scheme analysis'
        fig = plt.figure(title, figsize = (14,12))
        fig.suptitle(title, fontsize = 14 , fontweight = 'bold')
        
        ax = plt.subplot(221)
        xmin = min(kerogen.Ea) - 1
        xmax = max(kerogen.Ea) + 1
        ymin = 0
        ymax = max(kerogen.xi) + 50
        plt.axis([xmin, xmax, ymin, ymax])

        while T < maxT:
            for e in range(len(kerogen.Ea)):
                kerogen.xi[e] *= Arrhenius(kerogen.A , kerogen.Ea[e] , T , dts)

            # Compute HI
            HI = 0
            for c in kerogen.xi:
                HI += kerogen.dE * c

            for i in range(vitrinitetble.shape[0]):
                vitrinitetble[i][1] *= Arrhenius(1.0e13 , vitrinitetble[i][0] , T , dts)
            Tr = (85.0 - np.sum(vitrinitetble , axis = 0)[1]) / 100.0
            ezRo =  math.exp(-1.6 + 3.7 * Tr) / 100
            T2 = T - 15. * math.log(alpha / 2. , 10)

                # Initial heating : 5min @ 300C
            TR =  (1 - HI / kerogen.HI0)
            if TR < 1:
                RE_table = [kerogen.xi[e] for e in range(len(kerogen.Ea))]

                for e in range(len(kerogen.Ea)):
                    RE_table[e] *= Arrhenius(kerogen.A  , kerogen.Ea[e] , 300. , 300.)

                REYield = []
                for t in range(305,605,5):
                    for e in range(len(kerogen.Ea)):
                        RE_table[e] *= Arrhenius(kerogen.A  , kerogen.Ea[e] , t , 12)

                    REYielde = 0
                    for r in RE_table:
                        REYielde += kerogen.dE * r
                    REYield.append([t , REYielde])

                maxY = 0
                TMax1 = 0
                TMax2 = 0
                TMax3 = 0
                REYield = np.asarray(REYield)

                for k in range(1 , REYield.shape[0] - 1):
                    maxY = REYield[k - 1][1] - REYield[k][1]
                    if maxY > (REYield[k][1] - REYield[k+1][1]):
                        TMax1 = REYield[k - 1][0]
                        TMax2 = REYield[k + 0][0]
                        TMax3 = REYield[k + 1][0]
                        if k > 1 :
                            dY1 = REYield[k - 2][1] - REYield[k - 1][1]
                        else:
                            dY1 = 0
                        dY2 = REYield[k - 1][1] - REYield[k + 0][1]
                        dY3 = REYield[k + 0][1] - REYield[k + 1][1]
                        break

                matrix = [[TMax1 * TMax1 , TMax1 , 1],
                        [TMax2 * TMax2 , TMax2 , 1],
                        [TMax3 * TMax3 , TMax3 , 1]]
                matrix = np.asarray(matrix)

                try:
                    invmat = np.linalg.inv(matrix)
                    A = invmat[0][0] * dY1 + invmat[0][1] * dY2 + invmat[0][2] * dY3
                    B = invmat[1][0] * dY1 + invmat[1][1] * dY2 + invmat[1][2] * dY3
                    TMax = -B / (2 * A)
                    TMax = (TMax + 14) / 1.114 # Espitalie 1986
                except:
                    T += alpha * dt
                    continue

            result.append([T , HI , TR , ezRo , TMax , T2 , age])

            plt.xlabel('Activation Energy (kcal/mol)')
            plt.ylabel('Partial potential (mgHC/gC)')
            plt.title('Ea distribution vs. Maturity')
            if T in temptable:
                # label = str('% 0.2f' % float(100 * ezRo)) + '%Ro eq.'
                label = 'STS:' + str('% 0.0f' % float(T2)) + 'C'
                if kerogen.formalism == 'Pepper' :
                    ax.plot(kerogen.Ea, kerogen.xi, '-' , label = label)
                if kerogen.formalism == 'IFP' :
                    ax.bar(kerogen.Ea , kerogen.xi , label = label)
                ax.legend(fontsize = 8)

            T += alpha * dt
            age += dt
        ##### PLOT RESULTS #####
        result = np.asarray(result)
        result_df = pd.DataFrame(result, columns = labels)

        ax = plt.subplot(222)
        plt.xlabel('TMax (C)')
        plt.ylabel('Hydrogen Index (mgHC/gC)')
        plt.title('HI vs. TMax')
        xmin = 390
        xmax = 510
        ymin = 0
        ymax = 800
        plt.axis([xmin, xmax, ymin, ymax])
        a = result[:,1]
        b = result[:,4]

        ax.plot(b, a, '-' , b , a, lw = 2, c='k')
        # ax.fill(x, y, zorder=10)
        ax.grid(True, zorder=5)
        if user_data is not None:
            try:
                ax.scatter(user_data['TMAX'], user_data['HI'], c='r', s=50, marker='o')
            except:
                pass

        ax = plt.subplot(223)
        plt.xlabel('Thermal Stress (C)')
        plt.ylabel('Transformation ratio (/)')
        plt.title('TR vs. Maturity')
        xmin = 30
        xmax = 180
        ymin = max(scalar_start - 0.1, 0.0)
        ymax = scalar_end + 0.1
        plt.axis([xmin, xmax, ymin, ymax])
        y_values = result[:,2]*(scalar_end - scalar_start) + scalar_start
        a = result[:,5]
        b = y_values
        ax.plot(a, b, '-' , a , b, lw = 2, c='k')
        # ax.fill(x, y, zorder=10)
        ax.grid(True, zorder=5)
        if user_data is not None:
            try:
                ax.scatter(user_data['STS'], user_data['TR'], c='r', s=50, marker='o')
            except:
                pass

        ax = plt.subplot(224)
        plt.xlabel('Thermal Stress (C)')
        plt.ylabel('TMax (C)')
        plt.title('TMax vs. Maturity')
        xmin = 30
        xmax = 180
        ymin = 390
        ymax = 510
        plt.axis([xmin, xmax, ymin, ymax])
        a = result[:,4]
        b = result[:,5]
        ax.plot(b, a, '-' , b , a, lw = 2, c='k')
        # ax.fill(x, y, zorder=10)
        ax.grid(True, zorder=5)
        if user_data is not None:
            try:
                ax.scatter(user_data['STS'], user_data['TMAX'], c='r', s=50, marker='o')
            except:
                pass
        
        top_fig = plt.figure(figsize=(14,10))
        ax = top_fig.add_subplot(111)
        if top_fig_select == 'TR vs. STS':
            ax.set_xlabel('Thermal Stress (C)')
            ax.set_ylabel('Transformation ratio (/)')
            ax.set_title('TR vs. Maturity')
            xmin = 30
            xmax = 180
            ymin = max(scalar_start - 0.1, 0.0)
            ymax = scalar_end + 0.1
            ax.axis([xmin, xmax, ymin, ymax])
            y_values = result[:,2]*(scalar_end - scalar_start) + scalar_start
            a = result[:,5]
            b = y_values
            ax.plot(a, b, '-' , a , b, lw = 2, c='k')
            # ax.fill(x, y, zorder=10)
            ax.grid(True, zorder=5)
            if user_data is not None:
                try:
                    ax.scatter(user_data['STS'], user_data['TR'], c='r', s=50, marker='o')
                except:
                    pass
        elif top_fig_select == 'TMAX vs. STS':
            ax.set_xlabel('Thermal Stress (C)')
            ax.set_ylabel('TMAX (C)')
            ax.set_title('TMAX vs. STS')
            xmin = 30
            xmax = 180
            ymin = 390
            ymax = 510
            ax.axis([xmin, xmax, ymin, ymax])
            a = result[:,4]
            b = result[:,5]
            ax.plot(b, a, '-' , b , a, lw = 2, c='k')
            # ax.fill(x, y, zorder=10)
            ax.grid(True, zorder=5)
            if user_data is not None:
                try:
                    ax.scatter(user_data['STS'], user_data['TMAX'], c='r', s=50, marker='o')
                except:
                    pass
        elif top_fig_select == 'HI vs. TMAX':
            ax.set_xlabel('TMAX (C)')
            ax.set_ylabel('Hydrogen Index (mgHC/gC)')
            ax.set_title('HI vs. TMAX')
            xmin = 390
            xmax = 510
            ymin = 0
            ymax = 800
            ax.axis([xmin, xmax, ymin, ymax])
            a = result[:,1]
            b = result[:,4]

            ax.plot(b, a, '-' , b , a, lw = 2, c='k')
            # ax.fill(x, y, zorder=10)
            ax.grid(True, zorder=5)
            if user_data is not None:
                try:
                    ax.scatter(user_data['TMAX'], user_data['HI'], c='r', s=50, marker='o')
                except:
                    pass




    return result_df, fig, top_fig




class KerogenPepper:
    def __init__(self, name, A , Emean, s , HI0):
        self.name = name
        self.formalism = 'Pepper'
        self.A = A
        self.HI0 = None
        self.xi = None
        self.s = s
        self.Emean = Emean
        factor = 5.
        Ea = np.linspace(Emean - factor * s , Emean + factor * s , 250, retstep = True)
        self.Ea = Ea[0]
        self.dE = Ea[1]
        self.update_xi(HI0)

    def update_xi(self, HI0):
        self.HI0 = HI0
        self.xi = list(map(lambda e : HI0 *  ((1 / (self.s * math.sqrt(2 * math.pi))) * math.exp(-(math.pow(e - self.Emean , 2)) / (2 * math.pow(self.s , 2)))) , self.Ea))


# def compute(HI, OF, alpha):
    
#     of_select = of_dict[OF]
#     kero  = KerogenPepper(of_select['Name'] ,  of_select['A'] ,  of_select['E'] ,  of_select['s'] , float(HI))

#     primarycracking(kero , float(alpha))

#     with open("a.png", "rb") as image_file:
#         encoded_string = base64.b64encode(image_file.read())

#     return [{"type": "image", "data": {"alt": "could not compute", "src": "data:image/png;base64, " + encoded_string.decode('ascii')}}]

def compute_cracking(OF, HI, alpha, user_data=None, top_fig_select='TR vs. STS', scalar_start=0, scalar_end=1):
    '''
    Compute primary cracking for a given kerogen and heating rate

    Parameters:
    - OF (str) : Organo facies label ['A', 'B', 'C', 'DE', 'F']
    - HI0 (int) : initial Hydrogen Index (mgHC/gC)
    - alpha (float ): heating rate(C/Ma)

    Returns:
    - Pandas Dataframe with the simulated values
    - Matplotlib figure object
    '''
    plt.close()
    # of_select = of_dict[OF]
    of_select = OF
    kero  = KerogenPepper(of_select['Name'] ,  of_select['A'] ,  of_select['E'] ,  of_select['s'] , float(HI))
    with _lock:
        result_df, fig, top_fig = primarycracking(kero , float(alpha), user_data, top_fig_select, scalar_start, scalar_end)

    return result_df, fig, top_fig


def st_ui():
    st.set_page_config(layout = "wide")
    user_data = None
    user_file = st.sidebar.file_uploader("Data upload (Excel format, columns: 'STS', 'TR', 'HI', 'TMAX', no specific order, no mandatory columns)")
    if user_file is not None:
        user_data = pd.read_excel(user_file)
        user_data = user_data.fillna(-999)

    top_fig_select = st.sidebar.selectbox('Top plot selection', ['TR vs. STS', 'HI vs. TMAX', 'TMAX vs. STS'])
    st.title("Primary Cracking simulation")
    OF = st.sidebar.selectbox('OrganoFacies selection', ('A', 'B', 'C', 'DE', 'F'))
    A = st.sidebar.text_input("Frequency factor A (1e13 s-1)", value = of_dict[OF]['A']/1e13)
    E = st.sidebar.text_input("Mean activation energy (kJ/mol)", value = round(of_dict[OF]['E'],2))
    s = st.sidebar.text_input("Standard deviation (kJ/mol)", value = of_dict[OF]['s'])
    scalar_start = st.sidebar.text_input("Reaction startpoint", value = 0.)
    scalar_end = st.sidebar.text_input("Reaction endpoint", value = 1.)
    new_of = {'Name': "new_of", 'A': 1e13*float(A), 'E': float(E), 's': float(s)}
    HI = st.sidebar.text_input("Initial Hydrogen Index", value = 600)
    alpha = st.sidebar.slider("Heating rate", 0.1, 10., 2.)


    result_df, fig, top_fig = compute_cracking(new_of, HI, alpha, user_data, top_fig_select, float(scalar_start), float(scalar_end))
    
    st.header("Primary cracking dashboard")
    with _lock:
        buf = BytesIO()
        top_fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption='Primary Cracking dashboard')
        buf = BytesIO()
        fig.savefig(buf, format="png", bbox_inches='tight', transparent = True)
        st.image(buf, use_column_width=False, caption='Primary Cracking dashboard')
    
    st.header("Simulated data")
    st.dataframe(result_df)



if __name__ == '__main__':
    st_ui()