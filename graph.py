
import matplotlib.pyplot as plt

with open('../time.txt', 'r') as f:
    data = f.read()

    # Usuwamy znaki końca linii i dzielimy napis na sekcje
    sections = data.strip().split('\n')

    # Konwertujemy każdą sekcję na listę liczb całkowitych
    group1 = list(map(int, sections[0].strip(';').split(';')))
    group2 = list(map(int, sections[1].strip(';').split(';')))
    group3 = list(map(int, sections[2].strip(';').split(';')))
    group4 = list(map(int, sections[3].strip(';').split(';')))

plt.title('Czas rozwiązania układów, od wielokości macierzy')

plt.plot(group1, group2, label='Jacobi')

plt.plot(group1, group3, label='Gauss-Seidel')

plt.plot(group1, group4, label='Faktoryzacja LU')

plt.xlabel("Wielkość macierzy N")
plt.ylabel("Czas rozwiązania układu (ms)")

plt.legend()
plt.savefig('czas_N.png')
plt.show()

