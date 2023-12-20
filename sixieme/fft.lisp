(defpackage #:fft
  (:use #:cl)
  (:export #:fft
           #:schonhage-strassen))


(in-package #:fft)


(defmacro while (condition &body body)
  `(loop while ,condition
         do (progn ,@body)))


(defun ext-gcd (a b)
  (let ((s 0) (old-s 1) (r b) (old-r a)
        (quotient) (bezout-t))
    (while (not (zerop r))
      (setq quotient (floor old-r r))
      (psetq old-r r r (- old-r (* quotient r))
             old-s s s (- old-s (* quotient s))))
    (if (zerop b) (setq bezout-t 0)
        (setq bezout-t (floor (- old-r (* old-s a)) b)))
    (list old-r old-s bezout-t)))


(defun append-zeros (poly &optional (twice?))
  (when (not (vectorp poly)) (return-from append-zeros))
  (let ((len (length poly)))
    (when (and (not twice?) (zerop (logand len (1- len))))
      (return-from append-zeros poly))
    (let ((req-len))
      (do ((pow 2 (ash pow 1))) ((>= pow len) (setq req-len pow)))
      (when twice? (setq req-len (ash req-len 1)))
      (concatenate 'vector poly (make-array (- req-len len)
                                            :initial-element 0)))))


(defun rotate (poly len)
  (when (or (not (vectorp poly)) (not (zerop (logand len (1- len)))))
    (return-from rotate))
  (let ((oldest -1) (bit-len (round (log len 2)))
        (rev (make-array len :initial-element 0)))
    (do ((mask 1 (1+ mask))) ((= len mask) poly)
      (when (zerop (logand mask (1- mask)))
        (setq oldest (1+ oldest)))
      (setf (aref rev mask) (logior (aref rev (logxor mask (ash 1 oldest)))
                                    (ash 1 (- bit-len oldest 1))))
      (when (> (aref rev mask) mask)
        (rotatef (aref poly mask) (aref poly (aref rev mask)))))))


(defun fft (poly &optional mode)
  (when (not (vectorp poly)) (return-from fft))
  (let (len root roots alpha beta (divisor 1) (omega 4) m group)
    (setq poly (append-zeros poly)
          len (length poly)
          poly (rotate poly len)
          m (1+ (expt omega (ash len -1)))

          group (make-array len :initial-contents (loop for j from 0 below len
                                                        collect (mod (expt omega j) m))))
    (when mode
      (setq group (map 'vector #'(lambda (omega-pow)
                                   (mod (cadr (ext-gcd omega-pow m)) m)) group)))
    (format t "~% Значение модуля m: ~a.~%" m)
    (format t "~% Группа, порождаемая корнем омега ~a-ой степени: ~a.~%" len group)
    (do ((j 0 (1+ j))) ((= (round (log len 2)) j))
      (setq root (expt omega (ash (ash len -1) (- j))))
      (when (equal 'INV mode) (setq root (mod (cadr (ext-gcd root m)) m)))
      (setq roots (map 'vector #'(lambda (pow) (mod (expt root pow) m))
                       (loop for pow from 0 below (ash 1 j) collect pow)))
      (format t "~% Корни, полученные на текущей (~a) итерации: ~a.~%" (1+ j) roots)
      (do ((i 0 (+ i (ash 1 (1+ j))))) ((= len i))
        (do ((s 0 (1+ s))) ((= (ash 1 j) s))
          (setq alpha (aref poly (+ i s))
                beta (* (aref poly (+ i s (ash 1 j))) (aref roots s)))
          (setf (aref poly (+ i s)) (mod (+ alpha beta) m)
                (aref poly (+ i s (ash 1 j))) (mod (- alpha beta) m)))))
    (when (equal 'INV mode) (setq divisor (cadr (ext-gcd len m))))
    (map 'vector #'(lambda (cf) (mod (* cf divisor) m)) poly)))


(defun fft-handler (&optional mode)
  (format t "~% Введите коэффициенты многочлена от младшего к старшему через пробел (для ввода комплексного числа: #c(i j)):~2% Вводите: ")
  (let (input clear-input len-input cur res)
    (setq input (uiop:split-string (read-line) :separator " ")
          len-input (length input))
    (do ((j 0 (1+ j))) ((= j len-input) (setq clear-input (reverse clear-input)))
      (setq cur (nth j input))
      (handler-case (setq clear-input (cons (parse-integer cur) clear-input))
        (error ()
          (when (find #\( cur)
            (setq cur (concatenate 'string cur " " (nth (1+ j) input))
                  j (1+ j)
                  clear-input (cons (read-from-string cur) clear-input))))))
    (setq res (fft (make-array (length clear-input) :initial-contents clear-input) mode))
    (if mode
        (format t "~% Результат вычисления обратного быстрого преобразования Фурье для заданного многочлена:~2% ~a~%" res)
        (format t "~% Результат вычисления быстрого преобразования Фурье для заданного многочлена:~2% ~a~%" res))))


(defun schonhage-strassen-machinerie (&optional f-poly s-poly)
  (let* ((old-f-len (length f-poly)) (old-s-len (length s-poly))
         (sum-lens (+ old-f-len old-s-len)) (f-len old-f-len)
         (s-len old-s-len) res fft-f fft-s m
         (f-num (reduce #'+ (map 'vector #'(lambda (mult pow) (* mult (expt 10 pow)))
                                 f-poly (loop for pow from 0 below old-f-len collect pow))))
         (s-num (reduce #'+ (map 'vector #'(lambda (mult pow) (* mult (expt 10 pow)))
                                 s-poly (loop for pow from 0 below old-s-len collect pow)))))
    (if (> f-len s-len)
        (setq f-poly (append-zeros f-poly t) f-len (length f-poly)
              s-poly (concatenate 'vector s-poly
                                  (make-array (- f-len s-len) :initial-element 0)))
        (setq s-poly (append-zeros s-poly t) s-len (length s-poly)
              f-poly (concatenate 'vector f-poly
                                  (make-array (- s-len f-len) :initial-element 0))))
    (setq m (1+ (ash 1 (length f-poly))))
    (setq res (fft (map 'vector #'(lambda (f s) (mod (* f s) m))
                        (setq fft-f (fft f-poly)) (setq fft-s (fft s-poly))) 'INV))
    (format t "~% Результат применения быстрого преобразования Фурье к числу ~a, представленному в виде многочлена:~2% ~a.~%" f-num fft-f)
    (format t "~% Результат применения быстрого преобразования Фурье к числу ~a, представленному в виде многочлена:~2% ~a.~%" s-num fft-s)
    (reduce #'+ (map 'vector #'(lambda (cf pow)
                                 (* (expt 10 pow) (round (realpart cf))))
                     (subseq res 0 sum-lens)
                     (loop for i from 0 to (1- sum-lens) collect i)))))


(defun get-digits (num)
  (let ((digits (make-array 0 :adjustable t :fill-pointer 0)))
    (do () ((zerop num) digits)
      (vector-push-extend (mod num 10) digits)
      (setq num (floor num 10)))))


(defun schonhage-strassen ()
  (let (input res sign)
    (format t "~% Для умножения введите целые числа через пробел: ")
    (tagbody try-again
       (setq input (uiop:split-string (read-line) :separator " "))
       (handler-case
           (setq input (mapcar #'(lambda (num?) (parse-integer num?)) input)
                 input (cons 100 input)
                 sign (length (remove-if-not #'minusp input))
                 input (mapcar #'abs input))
         (error (err)
           (format t "~% В ходе выполнения программы была получена ошибка:~% ~a
 Введите числа для умножения снова: " err) (go try-again))))
    (setq res (reduce #'(lambda (f s)
                          (get-digits (schonhage-strassen-machinerie f s)))
                      (mapcar #'get-digits (subseq input 0 (1- (length input)))))
          res (schonhage-strassen-machinerie res (get-digits (car (last input))))
          res (/ (* res (setq sign (expt -1 (mod sign 2)))) 100))
    (format t "~% Результат произведения чисел равен: ~a.~%" res)
    (format t "~% Проверка: ~a = ~a: ~a.~%" res (setq input (apply #'* sign (cdr input))) (eql res input)) t))


(defun caller ()
  (let ((choice -1))
    (while t
      (format t "~% [1] -- Быстрое преобразование Фурье;
 [2] -- Обратное быстрое преобразование Фурье;
 [3] -- Алгоритм Шенхаге-Штрассена умножения натуральных чисел;
 [0] -- Выход.~%")
      (format t "~% Ваш выбор: ")
      (while (not (member (setq choice (read)) '(1 2 3 0)))
        (format t "~% Некорректный ввод! Попробуйте снова: "))
      (cond ((= 1 choice) (fft-handler))
            ((= 2 choice) (fft-handler 'INV))
            ((= 3 choice) (schonhage-strassen))
            (t (return-from caller t))))))


#|
(defun test ()
  (let ((times-incorrect 0) bad-number-pairs bad-numbers total res)
    (do ((i -100 (1+ i))) ((= i 100))
      (do ((j -100 (1+ j))) ((= j 100))
        (setq res (schonhage-strassen-machinerie (get-digits (abs i)) (get-digits (abs j))))
        (when (/= res (abs (* i j)))
          (incf times-incorrect)
          (setq bad-number-pairs (cons (list i j) bad-number-pairs)))))
    (setq bad-numbers (remove-duplicates (mapcar #'car bad-number-pairs))
          total (mapcar #'(lambda (bad-number) (list bad-number (length (remove-if-not #'(lambda (pair)
                                                                                           (= bad-number (car pair)))
                                                                                       bad-number-pairs))))
                        bad-numbers))
    (list total (apply #'+ (mapcar #'cadr total)))))
|#
