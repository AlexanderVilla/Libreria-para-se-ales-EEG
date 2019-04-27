import matplotlib.pyplot as pl
import numpy as np
from scipy import signal
import LinearFIR
from numpy import shape
import scipy 
from scipy.stats import kurtosis


def plot(senal, segundos, canal, nombre):
    Fs=250
    vector_tiempo = np.arange(0, len(senal)/Fs,1/Fs) 
    pl.figure(figsize=(10,5)) 
        #Se extrae canal 1 de la señal
    pl.grid('on')
    pl.title(nombre)
    pl.xlabel('Tiempo [s]')
    pl.ylabel('Amplitud[mV]') 
    pl.plot(vector_tiempo[0:Fs*segundos],senal[0:Fs*segundos, canal])
    pl.show()
    
def segmentacion(senal, segundos, Fs):
    epocas = senal.shape[1]//(Fs*segundos)
    print(epocas)
    senal = senal[:,0:epocas*Fs*segundos][:]
    segm = np.reshape(senal, (8 , Fs*segundos, epocas), order = 'F') # canales, muestras , epocas
    return(segm)

def valores_extremos(senal_segmentada, superior, inferior):
    canales, epocas = senal_segmentada.shape[0], senal_segmentada.shape[2] 
    senal_valores_extremos = senal_segmentada[:]
    epocas_malas = []
    for canal in range(canales):
        for epoca in range(epocas):
            maximo_aux = senal_segmentada[canal,:,epoca].max()
            minimo_aux = senal_segmentada[canal,:,epoca].min()
            if maximo_aux > superior or abs(minimo_aux) > abs(inferior):
                if epoca not in  epocas_malas:
                   epocas_malas.append(epoca)
    senal_valores_extremos = np.delete(senal_segmentada, epocas_malas, 2)
    print('Las épocas atípicas de la señal son: ', epocas_malas)
    
    #Para verificar
    maximo_new, minimo_new = senal_valores_extremos.max(), senal_valores_extremos.min()
    print(' Nueva Amplitud maxima: %d , Nueva amplitud mínima: %d' % (maximo_new, minimo_new))
    return(senal_valores_extremos)
    
def tendencia_lineal(data):
    lista_tendencias = []
    for canal in range(8):
        lista_tendencias.append(scipy.signal.detrend(data[:,canal]))
    return(np.array(lista_tendencias))

def improbabilidad(data):
    canales, epocas = data.shape[0], data.shape[2] 
    print(canales, epocas)
    Kurto= np.zeros((canales, epocas))
    for canal in range(canales):
        for epoca in range(epocas):
               Kurto[canal, epoca]=kurtosis(data[canal, : ,epoca])
    return(Kurto)

def patron_espectral(data, Fs, umbral):
    canales, epocas = data.shape[0], data.shape[2]
    epocas_malas = []
    for canal in range(canales):
        for epoca in range(epocas):
            f, Pxx = signal.welch(data[canal, :, epoca], Fs)
            nuevo_Pxx = Pxx- np.mean(Pxx)
            if  nuevo_Pxx.max()>umbral:
                if epoca not in  epocas_malas:
                   epocas_malas.append(epoca)
    print(epocas_malas)
    nuevo_arreglo = np.delete(data, epocas_malas, 2)
    return(nuevo_arreglo)
                
def welch(data):
    epocas = data.shape[2]
    lista = []
    for epoca in range(epocas):
            f, aux = signal.welch(data[:, :, epoca], Fs,'hamming')
            lista.append(aux )
    return(lista)
if __name__=="__main__":
    sec  = 20
    
    #============Cargar señales================================================
    senal = 'P1_RAWEEG_2018-11-15_FinProcedimiento_53min.txt'
    senalcolumnas = np.loadtxt(senal,delimiter=',',skiprows=6,usecols=[1,2,3,4,5,6,7,8])    
    Fs=250
    vector_tiempo = np.arange(0, len(senalcolumnas)/Fs,1/Fs) 
    
    plot(senalcolumnas, sec , 0, 'Señal Original')    
    #Parametros del filtro pasabajas
    N = 31# Orden del filtro (número de coeficientes)
    Wcbaja=1;
    Wcalta=50;
    Fny=Fs/2
    
    #Aplicacion del filtro 
    senalfiltrada = LinearFIR.eegfiltnew(senalcolumnas,Fs,Wcbaja,Wcalta)   
    plot(senalfiltrada, sec, 0, 'filtrada')
    
#    pl.figure(figsize=(10,5)) 
#    pl.grid('on')
#    pl.title('Señal Filtrada de las columnas')
#    pl.xlabel('Tiempo [s]')
#    pl.ylabel('Amplitud [V]')
#    pl.plot(vector_tiempo[0:10*Fs], senalfiltrada[0:10*Fs])
#    pl.show()
    
    #Aplicando nivel DC
    nivelDC=[0,100,200,300,400,500,600,700]
    senalfiltradaDC=senalfiltrada+nivelDC
    pl.figure(figsize=(10,5)) 
    pl.grid('on')
    pl.title('Señal Filtrada de las columnas con DC')
    pl.xlabel('Tiempo [s]')
    pl.ylabel('Amplitud [V]')
    pl.plot(vector_tiempo[0:sec*Fs],senalfiltradaDC[0:sec*Fs])
    pl.show()
    
     #==================Tendencia lineal=======================================
    senal_tendencia = tendencia_lineal (senalfiltrada)
    senal_tendencia = senal_tendencia.transpose()[:]
    plot(senal_tendencia, sec, 0, 'Tendencia lineal')
    
   # ===================Segmentacion===========================================  
    segundos = 2   
    senal =  senalfiltrada.transpose()
    senal_segmentada = segmentacion( senal , segundos, 250)
    
#    canal, muestras, epocas = senal_segmentada.shape[0],senal_segmentada.shape[1],senal_segmentada.shape[2]
#    senal_2 = np.reshape(senal_segmentada,(canal, muestras*epocas), order ='F')
##    nivelDC=[0,100,200,300,400,500,600,700]
##    senal_DC = senal_2.transpose()+nivelDC
##    vector_tiempo = np.arange(0, len(senal_DC)/Fs,1/Fs) 
##    pl.figure(figsize=(10,5)) 
##    pl.plot(vector_tiempo[0:sec*Fs], senal_DC[0:sec*Fs])
##    pl.title('')
##    pl.grid()
###
    
    #==================Valores Extremos========================================
    maximo, minimo = senal_segmentada.max(), senal_segmentada.min()
    print('Amplitud maxima: %d , Amplitud mínima: %d' % (maximo, minimo))
    
    superior, inferior= input("Ingrese un porcentaje de  umbral superior: " ), input("Ingrese un porcentaje de  umbral inferior: " )
    superior, inferior= (int(superior)/100)*maximo, (int(inferior)/100)*minimo
    senal_v_e = valores_extremos(senal_segmentada, superior, inferior)
    
    canal, muestras, epocas = senal_v_e.shape[0],senal_v_e.shape[1],senal_v_e.shape[2]
    senal_2 = np.reshape(senal_v_e,(canal, muestras*epocas), order ='F')
    
    nivelDC=[0,100,200,300,400,500,600,700]
    senal_DC = senal_2.transpose()+nivelDC
    vector_tiempo = np.arange(0, len(senal_DC)/Fs,1/Fs) 
    pl.figure(figsize=(10,5)) 
    pl.plot(vector_tiempo[0:sec*Fs], senal_DC[0:sec*Fs])
    pl.title('Valores Extremos')
    pl.xlabel('Tiempo [s]')
    pl.ylabel('Amplitud [mV]')
    pl.grid()

#    #======================Improbabilidad======================================      
    imp= improbabilidad(senal_v_e)
    nivelDCimp=[0,3,6,9,12,15,18,21]
    pl.figure(figsize=(10,5)) 
    pl.plot(imp.transpose()+nivelDCimp)
    pl.title('Curtosis')
    pl.grid()
    pl.legend(list(range(8)))
#    
#    #====================Patron espectral======================================
    pxx = welch(senal_v_e)
    pxx = np.array(pxx)
    pxx_max = pxx.max()
    pxx_porcentaje = input("Ingrese un porcentaje de  umbral positivo: " )
    umbral =(int(pxx_porcentaje)/100)*pxx_max
    senal_3 = np.reshape(senal_v_e,(canal, muestras*epocas), order ='F')
    espectro = patron_espectral(senal_3, Fs, umbral)
    pl.figure(figsize=(10,5))
    pl.plot(espectro)
    
    f, Pxx = signal.welch(senal_v_e , Fs,'hamming')
    nuevo_Pxx = Pxx- np.mean(Pxx)
    pl.figure(figsize=(10,5))
    pl.plot(Pxx[0,:,51])
    pl.plot(nuevo_Pxx)
    pl.grid()
