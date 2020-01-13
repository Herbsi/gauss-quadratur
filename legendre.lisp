(ql:quickload :iterate)
(ql:quickload :alexandria)
(use-package 'iter)
(ql:quickload :trivia)
(use-package 'trivia)
(in-optimizer :trivial)

(defun legendre (k)
  (declare
   (optimize (speed 3))
   (sb-ext:muffle-conditions cl:style-warning))
  "Returns a funktion that evaluates the k'th Legendre Polynomial at x"
  (match k
    (0 (lambda (x) 1))
    (1 (lambda (x) x))
    (k (lambda (x)
         (declare (double-float x))
         (iter
           (for l from 2 to k)
           (for x-2 previous x-1 initially 1) ; value of (l-2)th LP(x)
           (for x-1 previous x-0 initially x) ; value of (l-1)th LP(x)
           (for x-0 =                         ; value of lth LP(x)
                (/ (- (* (1- (* 2 l)) x x-1)
                      (* (1- l) x-2))
                   l))
           (finally (return x-0)))))))

(defun legendre-1st-deriv (k)
  "Returns a function that evaluates the first derivative
of the k'th Legendre Polynomial"
  (declare
   (optimize (speed 3))
   (sb-ext:muffle-conditions cl:style-warning))
  (match k
    (0 (lambda (x) 0))
    (1 (lambda (x) 1))
    (k (lambda (x)
         (declare (double-float x))
         (iter
           (for l from 2 to k)
           (for x-1 = (funcall (legendre (1- l)) x)) ; value of (l-1)th LP(x)
           (for x1-2 previous x1-1 initially 0)      ; value of (l-2)th LP'(x)
           (for x1-1 previous x1-0 initially 1)      ; value of (l-1)th LP'(x)
           (for x1-0 =                               ; value of l th LP'(x)
                (/ (- (* (1- (* 2 l)) (+ (* x x1-1) x-1))
                      (* (1- l) x1-2))
                   l))
           (finally (return x1-0)))))))

(defun fixed-point (f x0)
  (declare
   (double-float x0)
   (optimize (speed 3)))
  (iter
    (repeat 10000)
    (for y0 previous yk initially x0)
    (for yk = (funcall f y0))
    (until (< (abs (- yk y0)) 1.0d-15))
    (finally (return yk))))

(defun newton-method (f df x0)
  (flet ((newton-transform (g dg)
           (lambda (x)
             (- x (/ (funcall g x) (funcall dg x))))))
    (fixed-point (newton-transform f df) x0)))

(defun legendre-roots (n)
  (iter
    (for k from 0 to n)
    (for x0 = (cos (/ (* (+ (* 4 k) 3) pi)
                      (+ (* 4 n) 6))))
    (collect (newton-method (legendre (1+ n))
                            (legendre-1st-deriv (1+ n))
                            x0))))

(defun integration-weight (roots)
  (iter
    (for xk in roots)
    (with n = (1- (length roots)))
    (collect (/ (* 2 (- 1 (expt xk 2)))
                (* (expt (+ n 1) 2)
                   (expt (funcall (legendre n) xk) 2))))))

(defun gauss-quadratur (n)
  (let* ((roots (legendre-roots n))
         (weights (integration-weight roots)))
    (lambda (f)
      (iter
        (for root in roots)
        (for weight in weights)
        (sum (* (funcall f root) weight))))))

(defmacro main (number-of-supports)
  `(list ,@(mapcar
            (lambda (n)
              (let* ((int (funcall (gauss-quadratur n) (lambda (x) (log (+ x 2.0d0)))))
                     (err (abs (- int (- (log 27.0d0) 2)))))
                `(cons ,int ,err)))
            number-of-supports)))

;; TODO remove hard-coded lambda and actual value from main macro

(let ((result (time (main (2 4 8 16)))))
  (iter
    (for n in '(2 4 8 16))
    (for (int . err) in result)
    (format t "n: ~2,,d     int: ~,16E     err: ~,16E~%" n int err)))
