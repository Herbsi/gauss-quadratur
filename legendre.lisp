(ql:quickload :iterate)
(use-package :iter)
(ql:quickload :alexandria)
(use-package :alexandria)

(defun legendre (k)
  (declare
   (optimize (speed 3))
   (sb-ext:muffle-conditions cl:style-warning))
  "Returns a funktion that evaluates the k'th Legendre Polynomial at x"
  (case k
    (0 (lambda (x) 1))
    (1 (lambda (x) x))
    (otherwise
     (lambda (x)
       (declare (double-float x))
       (iter
        (for l from 2 to k)
        (for x-2 previous x-1 initially 1)          ; value of (l-2)th LP(x)
        (for x-1 previous x-0 initially x)          ; value of (l-1)th LP(x)
        (for x-0 =                                  ; value of lth LP(x)
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
  (case k
    (0 (lambda (x) 0))
    (1 (lambda (x) 1))
    (otherwise
     (lambda (x)
       (declare (double-float x))
       (iter
        (for l from 2 to k)
        (for x-1 = (funcall (legendre (1- l)) x))      ; value of (l-1)th LP(x)
        (for x1-2 previous x1-1 initially 0)           ; value of (l-2)th LP'(x)
        (for x1-1 previous x1-0 initially 1)           ; value of (l-1)th LP'(x)
        (for x1-0 =                                    ; value of l th LP'(x)
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

(defun integration-weights (roots)
  (iter
   (for xk in roots)
   (with n = (1- (length roots)))
   (collect (/ (* 2 (- 1 (expt xk 2)))
               (* (expt (+ n 1) 2)
                  (expt (funcall (legendre n) xk) 2))))))

(defmacro lookup (table n set-if-not-found)
  `(if (gethash ,n ,table)
       (gethash ,n ,table)
       (setf (gethash ,n ,table) ,set-if-not-found)))

(let ((cached-roots (make-hash-table))
      (cached-weights (make-hash-table)))
  (defun gauss-quadratur (n)
    (let* ((roots (lookup cached-roots n (legendre-roots n)))
           (weights (lookup cached-weights n (integration-weights roots))))
      (lambda (f)
        (iter
         (for root in roots)
         (for weight in weights)
         (sum (* (funcall f root) weight)))))))

(defmacro define-integrator (name fn int-value)
  `(defmacro ,name (supports)
     `(list ,@(mapcar
               (lambda (n)
                 (let* ((int (funcall (gauss-quadratur n) #',fn))
                        (err (abs (- int ,int-value))))
                   `(cons ,int ,err)))
               supports))))

(define-integrator herwig (lambda (x) (log (+ x 2.0d0))) (- (log 27.0d0) 2))
(define-integrator joe (lambda (x) (/ 1 (+ 2 x))) (log 3.0d0))

(defmacro main (int n-values)
  (alexandria:with-gensyms (gn gint gerr)
    `(iter
      (with result = (time (,int ,n-values)))
      (for ,gn in (list ,@n-values))
      (for (,gint . ,gerr) in result)
      (format t "n: ~2,,d     int: ~,16E     err: ~,16E~%" ,gn ,gint ,gerr))))

(main herwig (2 4 8 16 32))
(main joe (2 4 8 16 32))
