(defpackage #:lattice-methods
  (:use #:cl))


(in-package #:lattice-methods)


(defun gauss-lattice-reduction-handler ()
  (let (dim b-1 b-2)
    (format t "~% Введите размерность n вектора в базисе: ")
    (setq dim (read))
    (format t "~% Введите компоненты b_{1i} первого вектора базиса решётки через пробел: ")
    (setq b-1 (make-array dim :initial-contents
                          (mapcar #'parse-integer (uiop:split-string (read-line) :separator " "))))
    (format t "~% Введите компоненты b_{2i} второго вектора базиса решётки через пробел: ")
    (setq b-2 (make-array dim :initial-contents
                          (mapcar #'parse-integer (uiop:split-string (read-line) :separator " "))))
    (when (aux:is-linear-dependent? (make-array (list 2 dim) :initial-contents (list b-1 b-2)))
      (format t "~% Введённые векторы не являются базисом! Завершение алгоритма.~%")
      (return-from gauss-lattice-reduction-handler nil))
    (unless (<= (ops:norm-sqr b-1) (ops:norm-sqr b-2))
      (psetf b-1 b-2
             b-2 b-1)
      (format t "~% Векторы поменялись местами, т.к. должно выполняться ||b_1|| <= ||b_2||.~%"))
    (destructuring-bind (rb-1 rb-2) (lattice-reduction:gauss-lattice-reduction b-1 b-2)
      (format t "~% Редуцированный по Минковскому базис решётки имеет вид:~%
~4tb_1 = ~a,~%~4tb_2 = ~a.~%" rb-1 rb-2))))


(defun lll-handler ()
  (let (basis input reduced-basis)
    (format t "~% Введите размерность m вектора в решётке и размерность n базиса через пробел (m >= n): ")
    (destructuring-bind (dim-m dim-n) (mapcar #'parse-integer (uiop:split-string (read-line) :separator " "))
      (unless (>= dim-m dim-n)
        (format t "~% m должно быть больше или равно n! Завершение алгоритма.~%")
        (return-from lll-handler nil))
      (do ((j 0 (1+ j))) ((= dim-n j))
        (tagbody try-again
           (format t "~% Введите m компонент ~a-го вектора базиса через пробел: " (1+ j))
           (setq input (uiop:split-string (read-line) :separator " "))
           (when (/= dim-m (length input))
             (format t "~% Количество компонент вектора не равно m! Попробуйте снова.~%")
             (go try-again)))
        (setq basis (cons (make-array dim-m :initial-contents
                                      (mapcar #'parse-integer input)) basis)))
      (when (aux:is-linear-dependent? (make-array (list dim-n dim-m) :initial-contents basis))
        (format t "~% Введённые векторы не являются базисом! Завершение алгоритма.~%")
        (return-from lll-handler nil))
      (setq reduced-basis (lll:lll (make-array (list dim-n dim-m)
                                               :initial-contents (reverse basis))))
      (format t "~% LLL-приведённый базис решётки имеет вид:~% ")
      (loop for j from 0 below dim-n
            for slice = (aux:array-slice reduced-basis j)
            do (format t "~%~4tb_~a = ~a.~%" (1+ j) slice)))))



(defun linear-programming-handler ()
  (let (n input cs basis x-cfs sle xs cur-sol dot-res temp)
    (format t "~% Введите значение n (размерность решётки и размерность вектора в базисе): ")
    (setq n (read))
    (format t "~% Введите строки неравенств с учётом коэффициентов c_{1i} и c_{2i}.~%")
    (do ((j 0 (1+ j))) ((= n j))
      (tagbody try-again
         (format t "~% Введите коэффициенты ~a-ой строки неравенства через пробел: " (1+ j))
         (setq input (uiop:split-string (read-line)))
         (when (/= (+ 2 n) (length input))
           (format t "~% Некорректное количество коэффициентов в строке! Попробуйте ввести их снова:~%")
           (go try-again))
         (setq input (mapcar #'parse-integer input)
               cs (cons (list (first input) (car (last input))) cs)
               basis (cons (make-array n :initial-contents (subseq input 1 (1+ n))) basis)
               sle (cons input sle))))
    (setq cs (reverse cs)
          x-cfs (reverse basis)
          basis (lp::transpose (make-array (list n n) :initial-contents x-cfs)))
    (when (aux:is-linear-dependent? basis)
      (format t "~% Введённые векторы не являются базисом! Завершение алгоритма.~%")
      (return-from linear-programming-handler nil))
    (setq sle (make-array (list n (+ 2 n)) :initial-contents (reverse sle))
          xs (lp:linear-programming basis sle))
    (when (zerop (length xs))
      (format t "~% Решений не найдено.~%")
      (return-from linear-programming-handler))
    (format t "~% Найденные решения x_i = (x_{i1}, ..., x_{in}) имеют вид:~2%")
    (loop for j from 0 below (length xs)
          do (format t "~4tx_~a = ~a.~%" (1+ j) (nth j xs)))
    (format t "~% Проверка. Подставим полученные решения в исходные неравенства.~%")
    (do ((i 0 (1+ i))) ((= (length xs) i))
      (format t  "~% x_~a:" (1+ i))
      (setq cur-sol (nth i xs))
      (do ((j 0 (1+ j))) ((= n j) (terpri))
        (setq temp (nth j cs))
        (format t "~%~4t~a <= ~a <= ~a: ~a"
                (car temp) (setq dot-res (ops:dot cur-sol (nth j x-cfs)))
                (cadr temp) (<= (car temp) dot-res (cadr temp)))))))


(defun caller ()
  (let ((choice -1))
    (aux:while t
      (format t "~% [1] -- Алгоритм Гаусса редукции решёток размерности 2;
 [2] -- Алгоритм Ленстры-Ленстры-Ловаша (LLL-алгоритм);
 [3] -- Алгоритм решения задачи целочисленного программирования;
 [0] -- Выход.~%")
      (format t "~% Ваш выбор: ")
      (aux:while (not (member (setq choice (read)) '(1 2 3 0)))
        (format t "~% Некорректный ввод! Попробуйте снова: "))
      (cond ((= 1 choice) (gauss-lattice-reduction-handler))
            ((= 2 choice) (lll-handler))
            ((= 3 choice) (linear-programming-handler))
            (t (return-from caller t))))))
