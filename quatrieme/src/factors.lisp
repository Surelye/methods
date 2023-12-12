; Функции для компактного отображения списка множителей

(defun n-elts (elt n)
  (if (> n 1)
      (list n elt)
      elt))


(defun compr (elt n lst)
  (if (null lst)
      (list (n-elts elt n))
      (let ((next (car lst)))
        (if (eql next elt)
            (compr elt (1+ n) (cdr lst))
            (cons (n-elts elt n) (compr next 1 (cdr lst)))))))


(defun compress (x)
  (if (consp x)
      (compr (car x) 1 (cdr x))
      x))


; Алгоритм разложения числа n rho-методом Полларда

(defun rho-pollard-machinerie (n x-0 &optional (c 1))
  (when (aux::miller-rabin n) (return-from rho-pollard-machinerie 'PRIME))
  (let ((mapping (lambda (x) (mod (+ c (* x x)) n)))
        (a x-0) (b x-0) (q))
    (tagbody map
       (setq a (funcall mapping a)
             b (funcall mapping (funcall mapping b))
             q (gcd (- a b) n))
       (cond ((and (< 1 q) (< q n)) (return-from rho-pollard-machinerie
                                      (list q (aux::miller-rabin q))))
             ((= n q) (return-from rho-pollard-machinerie))
             (t (go map))))))


(defun rho-pollard-wrapper (n x-0)
  (let ((c 1) (head) (factor) (factors))
    (macro::while (zerop (logand n 1))
      (setq factors (cons 2 factors) n (ash n -1)))
    (setq x-0 (mod x-0 n))
    (macro::while (/= 1 n)
      (setq factor (rho-pollard-machinerie n x-0 c))
      (cond ((eql 'PRIME factor) (setq factors (cons n factors) n 1))
            ((cadr factor) (setq factors (cons (setq head (car factor)) factors)
                                 n (/ n head)))
            ((null factor) (macro::while (= (- n 2)
                                            (setq c (1+ (random (1- n)))))))
            (t (setq n (/ n (setq head (car factor)))
                     factors (append factors
                                     (rho-pollard-wrapper head (random head)))))))
    factors))


(defun rho-pollard (n x-0)
  (compress (sort (rho-pollard-wrapper n x-0) #'<)))


; Алгоритм разложения числа n (p-1)-методом Полларда

(defun primep (n)
  (when (= 2 n) (return-from primep t))
  (loop for fac from 2 to (isqrt n)
        never (zerop (mod n fac))))


(defun generate-factor-base (n)
  (when (< n 2) (return-from generate-factor-base (list 2)))
  (let ((base))
    (when (> n 2) (setq base (cons 2 base)))
    (do ((prime? 3 (+ 2 prime?))) ((> prime? n) (reverse base))
      (when (primep prime?) (setq base (cons prime? base))))))


(defun p-1-pollard-machinerie (n)
  (when (= 3 n) (return-from p-1-pollard-machinerie 3))
  (let* ((log-n (log n)) (base-bound (sqrt (exp (sqrt (* log-n (log log-n))))))
         (factor-base (generate-factor-base (floor base-bound)))
         (a (+ 2 (random (- n 3)))) (b (gcd a n)) (l) (ret-val))
    (when (> b 1) (return-from p-1-pollard-machinerie b))
    (dolist (prime factor-base)
      (setq l (floor (/ log-n (log prime)))
            a (aux::mod-expt a (aux::mod-expt prime l n) n)))
    (setq b (gcd (1- a) n))
    (when (/= n b) (setq ret-val b)) ret-val))


(defun p-1-pollard-wrapper (n)
  (let ((factor) (factors))
    (macro::while (zerop (logand n 1))
      (setq factors (cons 2 factors)
            n (ash n -1)))
    (macro::while (/= 1 n)
      (tagbody try-again
         (setq factor (p-1-pollard-machinerie n))
         (when (or (null factor) (= 1 factor))
           (if (aux::miller-rabin n)
               (progn (setq factors (cons n factors))
                      (return-from p-1-pollard-wrapper factors))
               (go try-again)))
         (if (aux::miller-rabin factor)
             (setq factors (cons factor factors))
             (setq factors (append (p-1-pollard-wrapper factor) factors)))
         (setq n (/ n factor)))) factors))


(defun p-1-pollard (n)
  (compress (sort (p-1-pollard-wrapper n) #'<)))


; Алгоритм разложения числа n методом непрерывных дробей

(defun make-factor-base (factors)
  (let (factor-base)
    (setq factors (reduce #'append factors))
    (dolist (factor factors factor-base)
      (if (listp factor)
          (when (evenp (car factor))
            (setq factor-base (adjoin (cadr factor) factor-base)))
          (when (and (not (member factor factor-base))
                     (> (count factor factors) 1))
            (setq factor-base (cons factor factor-base)))))))


(defun cont-fracs-machinerie (n)
  (let* ((sqrt-n (sqrt n)) (P-prev 1) (P-cur (floor sqrt-n)) (a-k (floor sqrt-n))
         (x-k (- sqrt-n a-k)) (P_k^2 (mod (* P-cur P-cur) n)) recip-x-k
         P_k^2-factors factor-base (iter 1) (k (floor (log n 2))) (half-n (ash n -1))
         vecs len-vecs cur-vec ks Ps s gammas t-val q)
    (flet ((factor (num)
             (when (> num half-n) (setq num (- num n)))
             (let* ((absed (abs num)) (factors (rho-pollard absed (random absed))))
               (if (minusp num)
                   (setq factors (cons -1 factors))
                   factors)))
           (make-vec (factors)
             (let ((vec (loop for j from 0 below (length factor-base) collect 0)))
               (dolist (factor factors vec)
                 (if (atom factor)
                     (setf (nth (position factor factor-base) vec) 1)
                     (setf (nth (position (cadr factor) factor-base) vec) (car factor))))))
           (is-B-smooth (factors)
             (every #'(lambda (factor)
                        (if (atom factor)
                            (member factor factor-base)
                            (member (cadr factor) factor-base))) factors)))
      (setq P_k^2-factors (list (factor P_k^2)))
      (tagbody factor-residues
         (setq recip-x-k (/ x-k)
               a-k (floor recip-x-k) x-k (- recip-x-k a-k))
         (psetq P-prev P-cur
                P-cur (mod (+ (* a-k P-cur) P-prev) n)
                Ps (cons P-cur Ps))
         (setq P_k^2 (mod (* P-cur P-cur) n)
               P_k^2-factors (cons (factor P_k^2) P_k^2-factors))
         (when (< iter k) (setq iter (1+ iter))
               (go factor-residues))
         (setq factor-base (sort (make-factor-base P_k^2-factors) #'<)))
      (dolist (factors P_k^2-factors)
        (when (is-B-smooth factors)
          (setq vecs (cons (make-vec factors) vecs))))
      (setq len-vecs (length vecs))
      (do ((i 0 (1+ i))) ((= i (1- len-vecs)) ks)
        (setq cur-vec (nth i vecs))
        (do ((j (1+ i) (1+ j))) ((= j len-vecs))
          (when (every #'zerop (mapcar #'(lambda (f s) (mod (+ f s) 2))
                                       cur-vec (nth j vecs)))
            (setq ks (cons (list i j) ks)))))
      (handler-case (setq ks (nth (random (length ks)) ks))
        (error ()
          (return-from cont-fracs-machinerie (list NIL (rho-pollard n (random n))))))
      (setq Ps (reverse Ps)
            s (mod (apply #'* (mapcar #'(lambda (k)
                                          (nth k Ps)) ks)) n)
            vecs (mapcar #'(lambda (k) (nth k vecs)) ks)
            gammas (mapcar #'(lambda (f s)
                               (ash (+ f s) -1))
                           (cdr (car vecs)) (cdr (cadr vecs)))
            t-val (apply #'* (mapcar #'(lambda (factor pow) (expt factor pow))
                                     (cdr factor-base) gammas))
            q (list (gcd (+ s t-val) n) (gcd (- s t-val) n))))))


(defun cont-fracs-wrapper (n)
  (let ((bound-iters 1000) (current-iter 0) factors)
    (macro::while (/= 1 n)
      (when (= bound-iters current-iter)
        (return-from cont-fracs-wrapper (cons n factors)))
      (destructuring-bind (f s) (cont-fracs-machinerie n)
        (unless f
          (return-from cont-fracs-wrapper s))
        (if (and (= 1 f) (= 1 s))
            (incf current-iter)
            (progn (cond ((= 1 f) (setq factors (cons s factors) n (/ n s)))
                         ((= 1 s) (setq factors (cons f factors) n (/ n f)))
                         (t (setq factors (cons f (cons s factors)) n (/ n f s))))
                   (setq current-iter 0)))))
    factors))


(defun cont-fracs (n)
  (compress (sort (cont-fracs-wrapper n) #'<)))


(defun print-factorizations (factorizations n)
  (flet ((print-factorization (factorization)
           (let (current-factor (len-fact (length factorization)))
             (format t "~d = " n)
             (do ((j 0 (1+ j))) ((= j len-fact))
               (setq current-factor (nth j factorization))
               (if (atom current-factor)
                   (format t "~d" current-factor)
                   (format t "~d^~d" (cadr current-factor) (car current-factor)))
               (format t (if (= j (1- len-fact))
                             ".~%"
                             " * "))))))
    (destructuring-bind (f s t-val) factorizations
      (when f
        (format t "~%Результат разложения числа rho-методом Полларда: ")
        (print-factorization f))
      (when s
        (format t "~%Результат разложения числа p-1-методом Полларда: ")
        (print-factorization s))
      (when t-val
        (format t "~%Результат разложения числа методом непрерывных дробей: ")
        (print-factorization t-val)))) t)


(defun caller ()
  (let ((choice -1) n f-res s-res t-res)
    (flet ((try-cont-fracs ()
             (handler-case (cont-fracs n)
               (error ()
                 (rho-pollard n (random n))))))
      (macro::while t
        (format t "~%[1] -- Rho-метод Полларда;
[2] -- (p-1)-метод Полларда;
[3] -- Метод цепных дробей;
[4] -- Все методы;
[0] -- Выход.~%")
        (format t "~%Ваш выбор: ")
        (macro::while (not (member (setq choice (read)) '(1 2 3 4 0)))
          (format t "Некорректный ввод! Попробуйте снова: "))
        (when (zerop choice)
          (return-from caller t))
        (format t "~%Введите число, которое необходимо разложить на множители: ")
        (setq n (parse-integer (read-line)))
        (cond ((= 1 choice) (setq f-res (rho-pollard n (random n))))
              ((= 2 choice) (setq s-res (p-1-pollard n)))
              ((= 3 choice) (setq t-res (try-cont-fracs)))
              (t (setq f-res (rho-pollard n (random n))
                       s-res (p-1-pollard n)
                       t-res (try-cont-fracs))))
        (print-factorizations (list f-res s-res t-res) n)
        (setq f-res nil s-res nil t-res nil)))))
