import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Parameter-parameter dari soal
Q1 = 250  # Debit sungai (m³/detik)
Q2 = 125  # Debit buangan (m³/detik)
C1 = 20    # BOD sungai (mg/l)
C2 = 800   # BOD buangan (mg/l)
DO_sungai = 6  # DO sungai (mg/l)
DO_buangan = 8  # DO buangan (mg/l)
suhu_sungai = 31  # Suhu sungai (°C)
suhu_buangan = 22  # Suhu buangan (°C)
K2_20 = 0.23  # Koefisien reaerasi pada 20 °C (hari⁻¹)
K1_20 = 3.00  # Koefisien deoksigenasi pada 20 °C (hari⁻¹)
theta1 = 1.047  # Faktor koreksi suhu untuk reaerasi
theta2 = 1.024  # Faktor koreksi suhu untuk deoksigenasi
DO_baku_mutu = 5  # Baku mutu DO (mg/l)
DO_sat = 8.2  # DO jenuh pada suhu campuran (diasumsikan)

# 1. Menghitung BOD, DO, dan suhu campuran
CO = (Q1 * C1 + Q2 * C2) / (Q1 + Q2)
DO_campuran = (Q1 * DO_sungai + Q2 * DO_buangan) / (Q1 + Q2)
suhu_campuran = (Q1 * suhu_sungai + Q2 * suhu_buangan) / (Q1 + Q2)

# 2. Menghitung K1 dan K2 pada suhu campuran
K1 = K1_20 * theta2**(suhu_campuran - 20)
K2 = K2_20 * theta1**(suhu_campuran - 20)

# 3. Menghitung defisit oksigen awal (Da) dan kritis (Dc)
Da = DO_sat - DO_campuran
Dc = DO_sat - DO_baku_mutu

# 4. Mendefinisikan persamaan Streeter-Phelps untuk mencari tc
def persamaan_streeter_phelps(t, *data):
    K1, K2, La, Da, Dc = data
    return ((K1 * La) / (K2 - K1)) * (np.exp(-K1 * t) - np.exp(-K2 * t)) + Da * np.exp(-K1 * t) - Dc

# 5. Mencari BOD ultimate (La) terlebih dahulu
La = CO  # Menggunakan BOD campuran sebagai perkiraan awal La

# 6. Mencari waktu kritis (tc) dan perbaikan La berdasarkan tc
tc_awal = 0.1  # Tebakan awal untuk tc
data = (K1, K2, La, Da, Dc)

tc = fsolve(persamaan_streeter_phelps, tc_awal, args=data)[0]

# Setelah mendapatkan tc, hitung La kembali menggunakan nilai tc
La = (K1 * Dc) / (K2 * (np.exp(-K2 * tc) - np.exp(-K1 * tc))) + Da * np.exp(-K1 * tc)

# 7. Menghitung BOD maksimum yang diperbolehkan
BOD_maks = La / np.exp(K1 * tc)

# 8. Menghitung efisiensi pengolahan
C2_akhir = ((Q1 + Q2) * BOD_maks - Q1 * C1) / Q2
efisiensi = ((C2 - C2_akhir) / C2) * 100

# 9. Membuat profil pencemaran
def streeter_phelps(t):
    return ((K1 * La) / (K2 - K1)) * (np.exp(-K1 * t) - np.exp(-K2 * t)) + Da * np.exp(-K1 * t)

t = np.linspace(0, 1, 100)  # Menghasilkan 100 titik waktu antara 0 dan 1 hari
Dt = streeter_phelps(t)
DO = DO_sat - Dt

# 10. Visualisasi
plt.figure(figsize=(10, 6))
plt.plot(t, Dt, label='Defisit Oksigen (Dt)')
plt.plot(t, DO, label='DO')
plt.axhline(y=DO_sat, color='r', linestyle='--', label='DO Jenuh')
plt.axhline(y=DO_baku_mutu, color='g', linestyle='--', label='Baku Mutu DO')
plt.xlabel('Waktu (hari)')
plt.ylabel('Konsentrasi (mg/l)')
plt.title('Profil Pencemaran')
plt.legend()
plt.grid(True)
plt.show()

# 11. Menampilkan hasil
print("BOD maksimum yang boleh dibuang ke sungai:", BOD_maks, "mg/l")
print("Efisiensi pengolahan yang diperlukan:", efisiensi, "%")
