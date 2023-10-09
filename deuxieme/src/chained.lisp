; Алгоритм разложения чисел в цепную дробь

(defun div (dividend divider)
  (multiple-value-bind (quotient)
      (floor dividend divider) quotient))


(defun get-chained-fraction (numer denom)
  (let ((chained-fraction))
    (aux::while (not (zerop (mod numer denom)))
      (psetq chained-fraction (cons (div numer denom) chained-fraction)
             numer denom
             denom (mod numer denom)))
    (reverse (cons (div numer denom) chained-fraction))))


; Алгоритмы приложения цепных дробей
; Решение линейных диофантовых уравнений ax + by = c

(defun diophant-handler (a b c x y)
  (format t "Решение диофантова уравнения ~ax - ~ay = ~a имеет вид:~%" a b c)
  (format t "    x = (-1)^k * cQ_{k-1} + bt = ~a + ~at;~%" x b)
  (format t "    y = (-1)^k * cP_{k-1} + at = ~a + ~at.~%" y a))


(defun diophant-machinerie (a b c)
  (let* ((p-prev 0) (q-prev 1) (p-cur 1) (q-cur 0)
         (chained-fraction (get-chained-fraction a b))
         (len-chained (1- (length chained-fraction)))
         (q-i) (factor) (x) (y))
    (do ((i 0 (1+ i))) ((= len-chained i))
      (setq q-i (nth i chained-fraction))
      (psetq p-cur (+ (* q-i p-cur) p-prev) p-prev p-cur
             q-cur (+ (* q-i q-cur) q-prev) q-prev q-cur))
    (setq factor (expt -1 (logand (1+ len-chained) 1))
          x (* factor c q-cur)
          y (* factor c p-cur))
    (diophant-handler a b c x y)))

    
; Решение линейных сравнений ax = b (mod m)

(defun solve-congruence-machinerie (a b m)
  (let ((gcd-a-m (gcd a m)) (chained-fraction)
        (len-chained) (x))
    (unless (zerop (mod b gcd-a-m))
      (return-from solve-congruence-machinerie))
    (setq a (/ a gcd-a-m) b (/ b gcd-a-m) m (/ m gcd-a-m)
          chained-fraction (get-chained-fraction m a)
          len-chained (1- (length chained-fraction)))
    (do ((i 0 (1+ i)) (p-cur 1 (+ (* (nth i chained-fraction) p-cur) p-prev))
         (p-prev 0 p-cur))
        ((= len-chained i)
         (progn (setq x (mod (* (expt -1 (logand len-chained 1)) p-cur b) m))
                (mapcar #'(lambda (mult) (+ x (* mult m)))
                        (loop for mult from 0 to (1- gcd-a-m) collect mult)))))))


; Алгоритм вычисления символа Лежандра

(defun is-square-residue (a p)
  (= 1 (aux::mod-expt a (ash (1- p) -1) p)))


(defun compute-legendre-machinerie (a p)
  (cond ((zerop (mod a p)) 0)
        ((is-square-residue a p) 1)
        (t -1)))


(defun compute-legendre (a p)
  (when (and (integerp a) (aux::miller-rabin p))
    (compute-legendre-machinerie a p)))


; Алгоритм вычисления символа Якоби

(defun compute-jacobi-machinerie (a b)
  (let ((r 1) (t-val) (c))
    (when (< a 0) (setq a (- a))
      (when (= 3 (mod b 4)) (setq r (- r))))
    (tagbody eliminate-evenness
       (setq t-val 0)
       (aux::while (zerop (logand a 1))
         (setq t-val (1+ t-val) a (ash a -1)))
       (when (= 1 (logand t-val 1))
         (when (or (= 3 (mod b 8)) (= 5 (mod b 8)))
           (setq r (- r))))
       (when (and (= 3 (mod a 4)) (= 3 (mod b 4)))
         (setq r (- r)))
       (setq c a a (mod b c) b c)
       (if (not (zerop a))
           (go eliminate-evenness)
           (return-from compute-jacobi-machinerie r)))))


(defun compute-jacobi (a b)
  (when (and (integerp a) (integerp b) (= 1 (logand b 1)) (> b 1))
    (if (= 1 (gcd a b))
        (compute-jacobi-machinerie a b) 0)))


; Алгоритмы извлечения квадратного корня в Zp

(defun find-k-i (a-i q p)
  (do ((k 0 (1+ k))) ((= 1 (aux::mod-expt a-i (* (expt 2 k) q) p)) k)))


(defun get-inv (a p)
  (cadr (aux::ext-gcd a p)))


(defun seq-sqrt-Zp (a p)
  (when (/= 1 (compute-legendre a p)) (return-from seq-sqrt-Zp))
  (let ((b) (k-i -1) (k-is) (r-i) (m 0) (q (1- p))
        (a-prev a) (a-cur a) (pow))
    (aux::while (zerop (logand q 1)) (setq m (1+ m) q (ash q -1)))
    (aux::while (/= -1 (compute-legendre (setq b (random p)) p)))
    (aux::while (not (zerop k-i))
      (setq  k-i    (find-k-i a-cur q p) k-is (cons k-i k-is))
      (psetq a-cur  (mod (* a-cur (aux::mod-expt b (ash 1 (- m k-i)) p)) p)
             a-prev a-cur))
    (setq k-is (cdr k-is) r-i (aux::mod-expt a-prev (ash (1+ q) -1) p))
    (do ((i (length k-is) (1- i))) ((= 0 i) r-i)
      (setq pow (ash 1 (- m (car k-is) 1))
            r-i (mod (* r-i (get-inv (aux::mod-expt b pow p) p)) p)
            k-is (cdr k-is)))))


; Вероятностный алгоритм Чипполы

(defun complex-*-mod (f s sqr p)
  (let ((ffc (car f)) (fsc (cadr f))
        (sfc (car s)) (ssc (cadr s)))
    (mapcar #'(lambda (c) (mod c p))
            (list (+ (* ffc sfc) (* fsc ssc sqr))
                  (+ (* ffc ssc) (* fsc sfc))))))


(defun lst-mod-m (lst m)
  (mapcar #'(lambda (c) (mod c m)) lst))


(defun complex-^-mod (base power sqr m)
  (setq base (lst-mod-m base m))
  (do ((product '(1 0)))
      ((zerop power) (car product))
    (do () ((oddp power))
      (setq base (lst-mod-m (complex-*-mod base base sqr m) m)
            power (ash power -1)))
    (setq product (lst-mod-m (complex-*-mod product base sqr m) m)
          power (1- power))))


(defun cipolla (n p)
  (when (or (not (integerp n)) (not (integerp p)) (not (aux::miller-rabin p))
            (< p 3) (/= 1 (compute-legendre n p))) (return-from cipolla))
  (let ((a) (sqr))
    (aux::while (/= -1 (compute-legendre (- (* (setq a (random p)) a) n) p)))
    (setq sqr (- (* a a) n))
    (complex-^-mod (list a 1) (ash (1+ p) -1) sqr p)))
