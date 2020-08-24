import numpy as np
import pprint
import matplotlib.pyplot as plt
import sys
import cmath
from scipy import signal
import scipy.spatial.distance as sd

pp = pprint.PrettyPrinter(indent=4)
np.set_printoptions(threshold=sys.maxsize)

################################################
#           FUNCTIONS                          #
################################################


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y


def decode_BPSK(phi):


    if phi > 0 and phi < np.pi/2:
        return 1
    else:
        return -1

def decode_OOK(bit):
    rho,phi = cart2pol(bit.real,bit.imag)

    if rho > 0.25:
        return 1
    else:
        return 0


def decode_4ASK(bit):
    # using euclidean distance:
    distance_1 = abs(sd.euclidean(0.75,bit))
    distance_2 = abs(sd.euclidean(0.25,bit))
    distance_3 = abs(sd.euclidean(-0.5,bit))
    distance_4 = abs(sd.euclidean(-0.25, bit))

    distances = [distance_1,distance_2,distance_3,distance_4]

    if min(distances) == distance_1:
        return 0.75
    elif min(distances) == distance_2:
        return 0.25
    elif min(distances) == distance_3:
        return -0.75
    else:
        return -0.25

def decode_D8PSK(bit):
    distance_1 = sd.euclidean(np.pi/2, bit)
    distance_2 = sd.euclidean(np.pi/4, bit)
    distance_3 = sd.euclidean(0, bit)
    distance_4 = sd.euclidean(7*np.pi/4, bit)
    distance_5 = sd.euclidean(3*np.pi/2, bit)
    distance_6 = sd.euclidean(5*np.pi/4, bit)
    distance_7 = sd.euclidean(np.pi, bit)
    distance_8 = sd.euclidean(3*np.pi/4, bit)

    distances = [distance_1,distance_2,distance_3,distance_4,distance_5,distance_6,distance_7,distance_8]
    if min(distances) == distance_1:
        return np.pi/2
    elif min(distances) == distance_2:
        return np.pi/4
    elif min(distances) == distance_3:
        return 0
    elif min(distances) == distance_4:
        return 7*np.pi/4
    elif min(distances) == distance_5:
        return 3*np.pi/2
    elif min(distances) == distance_6:
        return 5*np.pi/4
    elif min(distances) == distance_7:
        return np.pi
    else:
        return 3*np.pi/4

#################################################
#       READ samples from file                  #
#################################################


def load(filename):
    x = np.fromfile(filename, np.uint8) - np.float32(
        127.5)  # by subtracting a float32 the resulting array will also be float32
    return 8e-3 * x.view(np.complex64)


samples = load("New_samples.dat")

BPSK = samples[164710:189663]
ASK_4 = samples[135847:159476]
D8PSK = samples[255249:277916]
QAM_16 = samples[340366:365665]
BPSK = samples[164710:189663]
OOK = samples[105017:130590]

data_rate_BPSK = 97400
data_rate_ASK = 102709
data_rate_OOK = 94905
data_rate_D8PSK = 107151
data_rate_QAM = 96098

phase_shift_ASK_D8PSK = 1.229
phase_shift_BPSK = 3.2
phase_shift_OOK = -0.23


Current_scheme = BPSK
Current_data_rate = data_rate_BPSK
phase = phase_shift_BPSK
###########################################
#    Generate sampling wave  and Sample  #
##########################################

pi=np.pi
#sine_wave = np.sin(2*np.pi*f * (x/fs)+phase)

t_scale = np.linspace(0, (1/2.4e6)*len(Current_scheme), len(Current_scheme))
sine_wave = np.cos(2*pi*Current_data_rate*t_scale + phase)
#sine_wave= signal.square (2* pi * freq * t_scale + phase , duty = 0.1)

Current_scheme = Current_scheme/max(Current_scheme)
sampled_Magnitude_data = np.zeros(shape=len(Current_scheme), dtype=complex)
sin_index = []
for i in range(len(Current_scheme)-1):
    if sine_wave[i] > sine_wave[i+1] and sine_wave[i] > sine_wave[i-1] :
        sample_point = Current_scheme[i]
        sampled_Magnitude_data[i] = sample_point

Initial_sampled_Magnitude_data = sampled_Magnitude_data
####################################################
#       Phase correction                            #
#####################################################
#
#
# #######################################################
# #       0/180 correction                              #
# ########################################################
#
sampled_Magnitude_data_phase = sampled_Magnitude_data[sampled_Magnitude_data != 0]
print(len(sampled_Magnitude_data_phase))
#
phase = np.zeros(shape=len(sampled_Magnitude_data_phase))
rho = np.zeros(shape=len(sampled_Magnitude_data_phase))

rho,phase = cart2pol(sampled_Magnitude_data_phase.real,sampled_Magnitude_data_phase.imag)
last_phase = 0
phases = []
values = []
for j in range(len(phase)-1):
    error = phase[j]-phase[j-1]
    if (error+last_phase) > np.pi/2 and (error+last_phase)< 3*np.pi/2:
        current_phase = 3*np.pi/2
        phases.append(current_phase+error)
        values.append(current_phase)
        last_phase = np.pi
    else:
        current_phase = np.pi/2
        phases.append(current_phase+error)
        values.append(current_phase)
        last_phase = np.pi/2

phases.append(np.pi/2)
values.append(0)
sampled_Magnitude_data_phase.real , sampled_Magnitude_data_phase.imag = pol2cart(rho,values)

####################################################
#               OOK                                #
###################################################
# phase = np.zeros(shape=len(sampled_Magnitude_data_phase))
# rho = np.zeros(shape=len(sampled_Magnitude_data_phase))
# rho,phase = cart2pol(sampled_Magnitude_data_phase.real,sampled_Magnitude_data_phase.imag)
# phases = []
# for j in range(len(phase)-1):
#     error = phase[j] - phase[j - 1]
#     if rho[j] > 0.2 and error < np.pi/2:
#         current_phase = np.pi/2
#         phases.append(current_phase + error )
#     else:
#         current_phase = phase[j]
#         phases.append(current_phase)
#
# phases.append(0)
# sampled_Magnitude_data_phase.real , sampled_Magnitude_data_phase.imag = pol2cart(rho,phases)
####################################################
#               BPSK                               #
###################################################
# phase = np.zeros(shape=len(sampled_Magnitude_data_phase))
# rho = np.zeros(shape=len(sampled_Magnitude_data_phase))
# rho,phase = cart2pol(sampled_Magnitude_data_phase.real,sampled_Magnitude_data_phase.imag)
# last_phase = 0
# phases = []
# values = []
# for j in range(len(phase)-1):
#     error = phase[j]-phase[j-1]
#     if (error) > np.pi/2 and last_phase == 0:
#         current_phase = np.pi
#         phases.append(current_phase + error)
#         values.append(current_phase)
#         last_phase = np.pi
#     else:
#         current_phase = 0
#         phases.append(current_phase + error )
#         values.append(current_phase)
#         last_phase = 0
# phases.append(0)
# sampled_Magnitude_data_phase.real , sampled_Magnitude_data_phase.imag = pol2cart(rho,phases)
####################################################
#       Phase correction D8PSK                      #
####################################################
# phase_D8PSK = np.zeros(shape=len(sampled_Magnitude_data_phase))
# rho_D8PSK = np.zeros(shape=len(sampled_Magnitude_data_phase))
#
# rho_D8PSK, phase_D8PSK = cart2pol(sampled_Magnitude_data_phase.real,sampled_Magnitude_data_phase.imag)
#
# last_phase = 0
# phases_D8PSK = []
# values = []
# for j in range(len(phase_D8PSK)-1):
#     error = phase_D8PSK[j]-phase_D8PSK[j-1]
#     if (error+last_phase) >= np.pi/4 and (error+last_phase)< np.pi/2:
#         current_phase = last_phase + np.pi/4
#         phases_D8PSK.append(current_phase+error)
#         values.append(current_phase)
#     else:
#         current_phase = last_phase
#         phases_D8PSK.append(current_phase+error)
#
#     elif (error+last_phase) >= np.pi/6 and (error+last_phase)< np.pi/3 :
#         current_phase = np.pi/4
#         phases_D8PSK.append(current_phase+error)
#         last_phase = np.pi/4
#         values.append(current_phase)
#
#     elif (error+last_phase) >= np.pi/3 and (error+last_phase) <= 2*np.pi/3:
#         current_phase = np.pi / 2
#         phases_D8PSK.append(current_phase + error)
#         last_phase = np.pi / 2
#         values.append(current_phase)
#
#
#     elif (error+last_phase) > 2*np.pi/3 and (error+last_phase) <= 5*np.pi/6:
#         current_phase = 3*np.pi / 4
#         phases_D8PSK.append(current_phase + error)
#         last_phase = 3*np.pi/4
#         values.append(current_phase)
#
#
#
#     elif (error+last_phase) > 5*np.pi/6 and (error+last_phase) < 7*np.pi/6:
#         current_phase = np.pi
#         phases_D8PSK.append(current_phase + error)
#         last_phase = np.pi
#         values.append(current_phase)
#
#
#     elif (error+last_phase) >= 7*np.pi/6 and (error+last_phase) <= 4*np.pi/3:
#         current_phase = 5*np.pi / 4
#         phases_D8PSK.append(current_phase + error)
#         last_phase = 5*np.pi / 4
#         values.append(current_phase)
#
#
#     elif (error + last_phase) > 4 * np.pi / 3 and (error + last_phase) < 5 * np.pi / 3:
#         current_phase = 3 * np.pi / 2
#         phases_D8PSK.append(current_phase + error)
#         last_phase = 3 * np.pi / 2
#         values.append(current_phase)
#
#
#     else:
#         current_phase = 7 * np.pi / 4
#         phases_D8PSK.append(current_phase + error)
#         last_phase = 7 * np.pi / 4
#         values.append(current_phase)

# values.append(0)
# phases_D8PSK.append(0)
# sampled_Magnitude_data_phase.real , sampled_Magnitude_data_phase.imag = pol2cart(rho_D8PSK, phases_D8PSK)
# print (values)
# print (np.pi/2)
#####################################################
#               16 - QAM                            #
#####################################################

#####################################################
#           DECODE values                           #
#####################################################
#
# def Decode_DBPSK(bit1,bit2):
#     if np.abs((np.angle(bit2)-np.angle(bit1))) >= 0.5*np.pi:
#         return 1
#     else:
#         return -1

# #
# # Decoded_DBPSK = np.zeros(shape=len(sampled_Magnitude_data), dtype=complex)
# # for j in range (len(sampled_Magnitude_data)-1):
# #     Decoded_DBPSK[j] = Decode_DBPSK(sampled_Magnitude_data[j],sampled_Magnitude_data[j-1])
# #
Decoded_4ASK = np.zeros(shape=len(sampled_Magnitude_data_phase), dtype=complex)
for j in range (len(sampled_Magnitude_data_phase)):
    Decoded_4ASK[j] = decode_4ASK((sampled_Magnitude_data_phase[j]))

Decoded_4ASK = Decoded_4ASK[3:]
print (len(Decoded_4ASK))

# Decoded_BPSK = np.zeros(shape=len(phase),dtype=complex)
# for j in range (len(phase)):
#     Decoded_BPSK[j] = decode_BPSK(phase[j])
#
# # Decoded_OOK = np.zeros(shape=len(sampled_Magnitude_data_phase),dtype=complex)
# # for j in range(len(sampled_Magnitude_data_phase)):
# #     Decoded_OOK[j] = decode_OOK(sampled_Magnitude_data_phase.real[j])
# # #
# # Decoded_8DPSK = np.zeros(shape=len(sampled_Magnitude_data_phase),dtype=complex)
# # for j in range (len(sampled_Magnitude_data_phase)):
# #     Decoded_8DPSK[j] = decode_D8PSK(sampled_Magnitude_data_phase[j])
#
#
# ##################################
# #        READING ASCII           #
# ##################################
# #       OOK                     #
# # Decoded_OOK = Decoded_OOK[4:]
# # message_OOK = []
# # print (len(Decoded_OOK))
# # for bit in Decoded_OOK.real:
# #     if bit == 1:
# #         message_OOK.append(1)
# #     else:
# #         message_OOK.append(0)
# #
# # characters = ''
# # byte = ''
# #
# # for bit in message_OOK:
# #     if len(byte) < 8:
# #         byte += str(bit)
# #
# #     else:
# #         characters += chr(int(byte, 2))
# #         byte = ''
# #
# # print("BPSK: ", characters)
#
# #################################
# #       BPSK                    #
# ################################
# Decoded_BPSK = Decoded_BPSK[4:]
# message_BPSK = []
#
# for bit in Decoded_BPSK.real:
#     if bit == 1:
#         message_BPSK.append(1)
#     else:
#         message_BPSK.append(0)
#
# characters = ''
# byte = ''
#
# for bit in message_BPSK:
#     if len(byte) < 8:
#         byte += str(bit)
#
#     else:
#         characters += chr(int(byte, 2))
#         byte = ''
#
# print("BPSK: ", characters)
# ##############################
#   DPSK                     #
##############################
# print(len(Decoded_DBPSK))
# Decoded_DBPSK = Decoded_DBPSK[Decoded_DBPSK != 0]
# print(len((Decoded_DBPSK)))
# Decoded_DBPSK = values
#
# message = []
#
# for bit in Decoded_DBPSK:
#     if bit == 0:
#         message.append(1)
#     else:
#         message.append(0)
#
# characters = ''
# byte = ''
#
# for bit in message:
#     if len(byte) < 8:
#         byte += str(bit)
#
#     else:
#         characters += chr(int(byte, 2))
#         byte = ''
#
# print(characters)
##############################
#   4ASK                     #
##############################
message = []

for bit in Decoded_4ASK.real:
    if bit == 0.75:
        message.append("10")
      #  message.append(0)
    elif bit == 0.25:
        message.append("11")
     #   message.append(1)
    elif bit == -0.25:
        message.append("01")
      #  message.append(1)
    else:
        message.append("00")
       # message.append(0)

characters = ''
byte = ''

for bit in message:
    if len(byte) < 8:
        byte += str(bit)

    else:
        characters += chr(int(byte, 2))
        print (byte)
        byte = ''

print("4 ASK: ", characters)
print(len(message))
#
# #######################################
# #       Decode 8DPSK                  #
# ######################################
# message_8DPSK = []
#
# for bit in Decoded_8DPSK:
#     if bit == np.pi/2:
#         message_8DPSK.append("011")
#     elif bit == np.pi/4:
#         message_8DPSK.append("001")
#     elif bit == 0:
#         message_8DPSK.append("000")
#     elif bit == 3*np.pi / 2:
#         message_8DPSK.append("101")
#     elif bit == 7*np/pi/4:
#         message_8DPSK.append("100")
#     elif bit == 5*np.pi/4:
#         message_8DPSK.append("111")
#     elif bit == np.pi:
#         message_8DPSK.append("110")
#     else:
#         message_8DPSK.append("010")
#
# characters = ''
# byte = ''
#
# for bit in message_8DPSK:
#     if len(byte) < 8:
#         byte += str(bit)
#
#     else:
#         characters += chr(int(byte, 2))
#         byte = ''
#
# print("8 DPSK: ", characters)
# print(len(message_8DPSK))
# #######################################
#       Plot DATA                     #
#######################################

plt.figure(1)
plt.plot(sampled_Magnitude_data.real, sampled_Magnitude_data.imag,'.')
plt.plot(sampled_Magnitude_data.real[0:9],sampled_Magnitude_data.imag[0:9],'X')
plt.xlabel('Real')
plt.ylabel('Imaginary')
plt.title('Constellation Diagram With Symbol Synchronisation')


plt.figure(2)
plt.plot(np.abs(Current_scheme), zorder=2)
plt.plot(abs(Initial_sampled_Magnitude_data), '.',zorder=3)
plt.plot(sine_wave, zorder=1)
plt.xlabel('Frequency(Hz)')
plt.ylabel('Magnitude (dB)')
plt.title('Absolute magnitude of data with sampling Sinusoid')

plt.figure(3)
plt.plot(np.unwrap(np.angle(Current_scheme)))

plt.figure(4)
plt.plot(sampled_Magnitude_data_phase.imag, sampled_Magnitude_data_phase.real, 'k.')
plt.plot(sampled_Magnitude_data_phase.imag[0:2], sampled_Magnitude_data_phase.real[0:2], '.')

plt.xlabel('Real')
plt.ylabel('Imaginary')
plt.title('Constellation Diagram with carrier synchronisation')


plt.figure(5)
fftspec = np.fft.fftshift(np.fft.fft(np.abs(Current_scheme)))
fftfreq = np.fft.fftshift(np.fft.fftfreq(len(Current_scheme), d=1/2.4e6))
plt.plot(fftfreq,10*np.log10(fftspec))
plt.xlabel('Frequency(Hz)')
plt.ylabel('Magnitude (dB)')
plt.title('Power Spectral Density')



plt.show()
