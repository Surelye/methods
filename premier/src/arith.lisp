; Корректное отображение списка

(defparameter *proper-formatting*
  "~{~#[~;x_~a~;x_~a и x_~a~:;~@{x_~a~#[~;, и ~:;, ~]~}~]~}")


; Алгоритм Евклида

(defun euclidean-algorithm-machinerie (f-num s-num)
  (when (> s-num f-num) (psetq f-num s-num
                               s-num f-num))
  (aux::while (not (zerop (mod f-num s-num)))
    (psetq f-num s-num
           s-num (mod f-num s-num)))
  s-num)


(defun euclidean-algorithm (f-num s-num)
  (cond ((not (and (integerp f-num) (integerp s-num))) nil)
        ((zerop f-num) s-num)
        ((zerop s-num) f-num)
        (t (euclidean-algorithm-machinerie (abs f-num) (abs s-num)))))


; Расширенный алгоритм Евклида

(defun div (dividend divider)
  (multiple-value-bind (quotient)
      (floor dividend divider) quotient))


(defun extended-euclidean-algorithm-machinerie (f-num s-num)
  (let ((prev-u 1) (u 0)
        (prev-v 0) (v 1)
        (quotient) (flag))
    (when (> s-num f-num) (psetq f-num s-num
                                 s-num f-num
                                 flag t))
    (aux::while (not (zerop (mod f-num s-num)))
      (psetq quotient (div f-num s-num)
             f-num s-num
             s-num (mod f-num s-num))
      (psetq prev-u u
             u (- prev-u (* quotient u))
             prev-v v
             v (- prev-v (* quotient v))))
    (if flag (list s-num v u) (list s-num u v))))


(defun extended-euclidean-algorithm (f-num s-num)
  (cond ((not (and (integerp f-num) (integerp s-num))) nil)
        ((and (zerop f-num) (zerop s-num)) (list 0 0 0))
        ((zerop f-num) (list s-num 0 1))
        ((zerop s-num) (list f-num 1 0))
        (t (extended-euclidean-algorithm-machinerie (abs f-num) (abs s-num)))))


; Бинарный алгоритм Евклида

(defun find-e (f-num s-num)
  (do* ((e 0 (+ 1 e))
        (buff (* s-num (expt 2 e)) (* 2 buff)))
       ((and (<= buff f-num) (< f-num (* 2 buff))) e)))


(defun lsbgcd (f-num s-num)
  (let ((e nil) (t-val nil) (buff nil))
    (when (> s-num f-num) (psetq f-num s-num
                                 s-num f-num))
    (aux::while (not (zerop s-num))
      (setq e (find-e f-num s-num)
            buff (* (expt 2 e) s-num)
            t-val (min (- (* 2 buff) f-num)
                       (- f-num buff)))
      (if (<= t-val s-num)
          (setq f-num s-num
                s-num t-val)
          (setq f-num t-val)))
    f-num))


(defun binary-euclidean (a b)
  (when (zerop a) (return-from binary-euclidean b))
  (when (or (zerop b) (= a b)) (return-from binary-euclidean a))
  (when (or (= 1 a) (= 1 b)) (return-from binary-euclidean 1))
  (when (and (zerop (logand a 1)) (zerop (logand b 1)))
    (return-from binary-euclidean
      (ash (binary-euclidean (ash a -1) (ash b -1)) 1)))
  (when (and (zerop (logand a 1)) (= 1 (logand b 1)))
    (return-from binary-euclidean (binary-euclidean (ash a -1) b)))
  (when (and (= 1 (logand a 1)) (zerop (logand b 1)))
    (return-from binary-euclidean (binary-euclidean a (ash b -1))))
  (when (and (= 1 (logand a 1)) (= 1 (logand b 1)))
    (if (> b a)
        (return-from binary-euclidean (binary-euclidean (ash (- b a) -1) a))
        (return-from binary-euclidean (binary-euclidean (ash (- a b) -1) b)))))


; Алгоритмы решения систем сравнений
; Китайская теорема об остатках при решении систем сравнений

(defun generate-new-congruences (M congruences)
  (let ((congruence-coeffs (mapcar #'(lambda (congruence)
                                       (/ M (cadr congruence)))
                                   congruences)))
    (mapcar #'cons congruence-coeffs congruences)))


(defun solve-congruence (coef res modulo)
  (let ((x (mod (cadr (extended-euclidean-algorithm coef modulo)) modulo)))
    (mod (* x res) modulo)))


(defun solve-congruence-system-machinerie (congruences)
  (let* ((M (apply #'* (mapcar #'cadr congruences)))
         (new-congruences (generate-new-congruences M congruences))
         (xs (mapcar #'(lambda (congruence)
                         (destructuring-bind (coef res modulo) congruence
                           (solve-congruence (mod coef modulo) res modulo)))
                     new-congruences)))
    (mod (apply #'+ (mapcar #'(lambda (congruence x) (* (car congruence) x))
                            new-congruences xs)) M)))


; Алгоритма Гарнера

(defun garners-algorithm-machinerie (congruences)
  (let* ((k (length congruences))
         (cs (loop for i from 2 to k collect 1))
         (coeffs (mapcar #'car congruences))
         (mods (mapcar #'cadr congruences))
         (u) (a-i) (m-i) (m-j) (c-i) (x) (i-pred))
    (do ((i 2 (1+ i))) ((> i k))
      (do ((j 1 (1+ j))) ((> j (1- i)))
        (setq i-pred (1- i) m-i (nth i-pred mods)
              m-j (nth (1- j) mods) c-i (nth (1- i-pred) cs)
              u (mod (cadr (extended-euclidean-algorithm m-j m-i)) m-i))
        (setf (nth (1- i-pred) cs) (mod (* u c-i) m-i))))
    (setq u (car coeffs)
          x u)
    (do ((i 2 (1+ i))) ((> i k) x)
      (setq i-pred (1- i) a-i (nth i-pred coeffs)
            c-i (nth (1- i-pred) cs) m-i (nth i-pred mods)
            u (mod (* (- a-i x) c-i) m-i)
            x (+ x (apply #'* u (subseq mods 0 i-pred)))))))


; Метод Гаусса решения систем линейных уравнений над конечными полями

(defun inv-mod (a p)
  (cadr (extended-euclidean-algorithm a p)))


(defun eliminate-a-n (sub-from to-sub n f-order)
  (let* ((el-coef (nth n sub-from))
         (to-sub-multed (mapcar #'(lambda (to-sub-elt)
                                    (mod (* to-sub-elt el-coef) f-order))
                                to-sub)))
    (mapcar #'(lambda (sub-from-elt to-sub-elt)
                (mod (- sub-from-elt to-sub-elt) f-order))
            sub-from to-sub-multed)))


(defun remove-trivials (lin-sys)
  (remove-if #'(lambda (equation)
                 (null (remove-if #'zerop equation)))
             lin-sys))


(defun find-nz-idx (coeffs &optional (fn #'not) (idx 0)) ; fn -> #'not or #'eval
  (cond ((null coeffs) nil)
        ((funcall fn (zerop (car coeffs))) idx)
        (t (find-nz-idx (cdr coeffs) fn (1+ idx)))))


(defun comp-routine (f s)
  (let* ((routine #'(lambda (lst) (find-nz-idx lst #'eval)))
         (f-res (funcall routine f)) (s-res (funcall routine s)))
    (cond ((null f-res) t)
          ((null s-res) nil)
          (t (> f-res s-res)))))


(defun gauss-method (lin-sys f-order)
  (let* ((lin-sys-copy (sort (remove-trivials lin-sys) #'comp-routine))
         (m (length lin-sys-copy)) (nz-idx) (nz-idx-eq) (a-is))
    (do ((i 0 (1+ i))) ((> i (1- m)))
      (setq a-is (mapcar #'(lambda (equation) (nth i equation))
                         (subseq lin-sys-copy i m))
            nz-idx (find-nz-idx a-is #'not))
      (when (null nz-idx)
        (format t "~%Решение не может быть найдено.~%")
        (return-from gauss-method))
      (setq nz-idx (+ i nz-idx))
      (setq nz-idx (+ i nz-idx)
            nz-idx-eq (nth nz-idx lin-sys-copy)
            nz-idx-eq (mapcar #'(lambda (a-i-j)
                                  (mod (* a-i-j (inv-mod (nth i nz-idx-eq)
                                                         f-order))
                                       f-order)) nz-idx-eq))
      (setf (nth nz-idx lin-sys-copy) nz-idx-eq)
      (do ((j 0 (1+ j))) ((> j (1- m)))
        (when (and (/= i j) (/= 0 (nth i (nth j lin-sys-copy))))
          (setf (nth j lin-sys-copy)
                (eliminate-a-n (nth j lin-sys-copy) nz-idx-eq i f-order))))
      (setq lin-sys-copy (remove-trivials lin-sys-copy)
            m (length lin-sys-copy)))
    lin-sys-copy))


; Обёртки


(defun parse-int (input-string &optional (start 0) (end (length input-string)))
  (handler-case (parse-integer input-string :start start :end end)
    (error (err)
      (format t "В ходе выполнения программы была получена ошибка:~%~a~%" err)
      (values nil))))


(defun construct-solution (lin-sys field-order to-input)
  (let ((binded (length lin-sys)) (input-string) (xs) (free-xs) (cur-eq))
    (format t "Введите значения свободных неизвестных (")
    (format t *proper-formatting*
            (loop for j from (1+ to-input) to (1- (length (car lin-sys))) collect j))
    (format t ") через пробел: ")
    (tagbody try-again
       (setq input-string (uiop:split-string (read-line) :separator " "))
       (handler-case (setq input-string (mapcar #'(lambda (num?)
                                                  (mod (parse-integer num?)
                                                       field-order))
                                              input-string))
         (error (err)
           (format t "В ходе выполнения программы была получена ошибка:~%~a
Введите значения свободных неизвестных снова: " err)
           (go try-again))))
    (do ((i 0 (1+ i))) ((= i binded) xs)
      (setq cur-eq (nth i lin-sys))
      (setq free-xs (mapcar #'(lambda (coef x-val) (* (- coef) x-val))
                            (subseq cur-eq binded) input-string)
            xs (cons (apply #'+ (car (last cur-eq)) free-xs) xs)))
    (format t "Вектор частного решения принимает вид: (~{~a~^, ~}).~%"
            (mapcar #'(lambda (x)
                        (mod x field-order))
                    (append (reverse xs) input-string)))))


(defun gauss-printer (lin-sys field-order)
  (let* ((solved (gauss-method lin-sys field-order)) (cur-eq)
         (num-eq (length solved)) (num-vars (1- (length (car lin-sys))))
         (to-input (- num-vars num-eq)))
    (when (null solved) (format t "~%") (return-from gauss-printer))
    (format t "~%Разрешённая матрица имеет вид:~%")
    (do ((i 0 (1+ i))) ((= i num-eq) (format t "~%"))
      (setq cur-eq (nth i solved))
      (do ((j 0 (1+ j))) ((= j num-vars) (format t " = ~a;~%" (car (last cur-eq))))
        (if (= j (1- num-vars)) (format t "~a * x_~a" (nth j cur-eq) (1+ j))
            (format t "~a * x_~a + " (nth j cur-eq) (1+ j)))))
    (if (= num-eq num-vars)
        (progn (format t "Таким образом, СЛУ имеет решение:~%")
               (do ((j 0 (1+ j))) ((= j num-vars) (format t "~%"))
                 (if (= j (1- num-vars))
                     (format t "    x_~a = ~a.~%" (1+ j) (car (last (nth j solved))))
                     (format t "    x_~a = ~a;~%" (1+ j) (car (last (nth j solved)))))))
        (progn (format t "Общее решение системы имеет вид:~%")
               (do ((i 0 (1+ i))) ((= i num-eq) (format t "~%"))
                 (format t "x_~a = " (1+ i))
                 (setq cur-eq (nth i solved))
                 (do ((j to-input (1+ j))) ((= j num-vars)
                                          (format t "~a;~%" (car (last (nth i solved)))))
                   (when (and (/= 0 (mod (- (nth j cur-eq)) field-order)) (/= (1+ i) (1+ j)))
                     (format t "~a * x_~a + " (mod (- (nth j cur-eq)) field-order) (1+ j)))))
               (format t "Хотите получить частное решение системы?")
               (when (y-or-n-p) (construct-solution solved field-order to-input))
               (format t "~%")))))


(defun gauss-handler ()
  (let ((num-eq) (num-vars) (field-order) (lin-sys) (input-line))
    (format t "~%Введите количество уравнений, входящих в СЛУ: ")
    (aux::while (not (and (integerp (setq num-eq (read))) (> num-eq 0)))
      (format t "Некорректное значение количества уравнений! Попробуйте снова: "))
    (format t "Введите количество переменных, входящих в каждое уравнение (без учёта свободного к/ф): ")
    (aux::while (not (and (integerp (setq num-vars (read))) (>= num-vars num-eq)))
      (format t "Некорректное значение количества переменных! Попробуйте снова: "))
    (format t "Введите порядок конечного поля: ")
    (aux::while (not (and (integerp (setq field-order (read))) (aux::miller-rabin field-order)))
      (format t "Некорректное значение порядка конечного поля! Попробуйте снова: "))
    (do ((m 1 (1+ m))) ((> m num-eq) (setq lin-sys (reverse lin-sys)))
      (format t "Введите коэффициенты a_~aj и свободный член строки #~a через пробел: " m m)
      (tagbody try-again
         (setq input-line (uiop:split-string (read-line) :separator " "))
         (handler-case (setq input-line (mapcar #'parse-integer input-line))
           (error (err)
             (format t "В ходе выполнения программы была получена ошибка:~%~a
Попробуйте ввести строку #~a снова: " err m) (go try-again)))
         (when (< (length input-line) (1+ num-vars))
           (format t "Количество введённых в строке #~a параметров недостаточно! Попробуйте ввести её снова: " m)
           (go try-again)))
      (setq lin-sys (cons (mapcar #'(lambda
                                        (num) (mod num field-order)) input-line)
                          lin-sys)))
    (gauss-printer lin-sys field-order)))


(defun sle-printer (congruences)
  (let ((m (apply #'* (mapcar #'cadr congruences)))
        (chinese (solve-congruence-system-machinerie congruences))
        (garners (garners-algorithm-machinerie congruences)))
    (format t "~%Результат решения системы линейных сравнений при помощи китайской теоремы об остатках:
    x = ~a -- наименьшее целочисленное решение. Остальные решения имеют вид: x + ~a * n, где n из Z.~%"
            chinese m)
    (format t "Результат решения системы линейных сравнений при помощи алгоритма Гарнера:
    x = ~a -- наименьшее целочисленное решение. Остальные решения имеют вид: x + ~a * n, где n из Z.~%~%"
            garners m)))


(defun sle-handler ()
  (let ((k) (congruences) (input-string) (pos) (a-j) (m-j)
        (m-js))
    (format t "~%Введите количество линейный сравнений: ")
    (aux::while (not (and (integerp (setq k (read))) (> k 0)))
      (format t "Некорректное количество линейных сравнений! Попробуйте снова: "))
    (do ((j 1 (1+ j))) ((> j k) (setq congruences (reverse congruences)))
      (format t "Введите значения a_~a и m_~a из соотношения x = a_~a (mod m_~a) через строку: "
              j j j j)
      (tagbody try-again
        (setq input-string (read-line)
               pos (position #\SPACE input-string))
         (when (null pos) (format t "Некорректный ввод! Попробуйте ввести строку #~a соотношения снова: " j)
               (go try-again))
         (aux::while (null (setq a-j (parse-int input-string 0 pos)))
           (format t "Некорректное значение параметра a_~a! Попробуйте ввести строку #~a соотношения снова: " j j)
           (go try-again))
         (aux::while (or (null (setq m-j (parse-int input-string (1+ pos))))
                           (remove-if #'(lambda (checked-m-j) (= 1 (gcd m-j checked-m-j))) m-js))
           (format t "Некорректное значение параметра m_~a! Попробуйте ввести строку #~a соотношения снова: " j j)
           (go try-again))
           (setq congruences (cons (list a-j m-j) congruences)
                 m-js (cons m-j m-js))))
    (sle-printer congruences)))


(defun euclidean-printer (a b)
  (let ((euc-res (euclidean-algorithm a b))
        (ext-euc-res (extended-euclidean-algorithm a b))
        (bin-euc-res (binary-euclidean a b)))
    (format t "Результат выполнения обычного алгоритма Евклида:
    gcd(~a, ~a) = ~a.~%" a b euc-res)
    (format t "Результат выполнения расширенного алгоритма Евклида:
    ~a * x + ~a * y = gcd(~a, ~a);
    gcd(~a, ~a) = ~a; x = ~a; y = ~a.~%"
            a b a b a b (car ext-euc-res) (cadr ext-euc-res) (caddr ext-euc-res))
    (format t "Результат выполнения бинарного алгоритма Евклида:
    gcd(~a, ~a) = ~a.~%~%" a b bin-euc-res)))


(defun euclidean-handler ()
  (let ((a) (b) (input-string) (pos))
    (format t "~%Введите значение коэффициентов a и b через пробел: ")
    (tagbody try-again
       (setq input-string (read-line)
             pos (position #\SPACE input-string))
       (when (null pos) (format t "Некорректный ввод! Попробуйте снова: ")
             (go try-again))
       (aux::while (null (setq a (parse-int input-string 0 pos)))
         (format t "Некорректное значение параметра a! Попробуйте снова: ")
         (go try-again))
       (aux::while (null (setq b (parse-int input-string (1+ pos))))
         (format t "Некорректное значение параметра b! Попробуйте снова: ")
         (go try-again)))
    (euclidean-printer a b)))


(defun caller ()
  (let ((choice -1))
    (aux::while t
      (format t "[1] -- Алгоритм Евклида (обычный, расширенный, бинарный);
[2] -- Алгоритмы решения сравнений (китайская теорема, метод Гарнера);
[3] -- Метод Гаусса решения СЛУ над конечными полями;
[0] -- Выход.~%")
      (format t "Ваш выбор: ")
      (aux::while (not (member (setq choice (read)) '(1 2 3 0)))
        (format t "Некорректный ввод! Попробуйте снова: "))
      (cond ((= 1 choice) (euclidean-handler))
            ((= 2 choice) (sle-handler))
            ((= 3 choice) (gauss-handler))
            (t (return-from caller))))))
