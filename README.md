# Конкурсное задание по суперкомпьютерной практике

## Формулировка задачи диффузии
Рассматривается решение нестационарного уравнения диффузии с неизвестной U.


## Аналитическое решение
<p align="center">
  <img src="data/images/analytical.png">
</p>

## Start code

```bash
pip install -r requirements.txt
mkdir build
cd build
cmake ..
make
./main 5 5 5
```