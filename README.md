# Программа расчёта положения центра инерции для молекул
Данная программа на C++ позволяет использовать метод молекулярной динамики с межмолекулярным взаимодейтсвием на основе потенциала Леннард-Джонса с упрощённым граничным условием (абсолютно упругое столкновение; нет термостата, периодических граничных условий и поршня), а также вычислять центр масс системы в данный момент времени. В программе для визуализации перемещений частиц используется модуль на Python, требующий numpy и matplotlib.
Пример визуализации полученных перемещений:
![image](https://github.com/alexdtat/CenterOfMass/assets/57017816/79d5d75a-6863-46f7-a089-b4a4e53c5597)
