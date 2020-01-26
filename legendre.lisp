(ql:quickload :herwigs-cl-utilities)
(ql:quickload :iterate)
(use-package :iter)
(ql:quickload :alexandria)

(defmacro make-legendre (name &key ((:nth-derivative n) 0))
  (alexandria:with-gensyms (x k l)
    (flet ((init-2 (k)
             (case k
               (0 1)
               (otherwise 0)))
           (init-1 (k)
             (case k
               (0 x)
               (1 1)
               (otherwise 0))))
      `(defun ,name (,k)
         (declare (sb-ext:muffle-conditions style-warning)
                  (optimize (speed 3)))
         (case ,k
           (0 (lambda (,x) ,(init-2 n)))
           (1 (lambda (,x) ,(init-1 n)))
           (otherwise
            (lambda (,x)
              (declare (double-float ,x))
              (iter (for ,l from 2 to ,k)
                ,@(let ((syms (herwig:group (herwig:map1-n #'gensym (* 3 (1+ n))) 3)))
                    (iter
                      (for k upfrom 0)
                      (for (x-0 x-1 x-2) in syms)
                      (appending
                       `((for ,x-2 previous ,x-1 initially ,(init-2 k))
                         (for ,x-1 previous ,x-0 initially ,(init-1 k))
                         (for ,x-0 = (/ (- (* (- (* 2.0d0 ,l) 1.0d0) (+ (* ,x ,x-1) ,@sofar))
                                           (* (- ,l 1.0d0) ,x-2)) ,l))) into acc)
                      (collect x-1 into sofar at start)
                      (finally (return (append acc
                                               `((finally (return ,(nth (* 3 n)
                                                                        (alexandria:flatten syms))))))))))))))))))

(make-legendre legendre-0)
(make-legendre legendre-1 :nth-derivative 1)

(defun fixed-point (f x0)
  (declare
   (double-float x0)
   (optimize (speed 3)))
  (flet ((good-enough? (xk xk-1)
           (< (abs (- xk xk-1)) 1.0d-15)))
    (iter
      (repeat 100)
      (for y0 previous yk initially x0)
      (format t "~& ~a ~%" y0)
      (for yk = (funcall f y0))
      (if (good-enough? yk y0) (return yk))
      (finally (return yk)))))

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
    (collect (newton-method (legendre-0 (1+ n))
                            (legendre-1 (1+ n))
                            x0))))

(defun integration-weights (roots)
  (iter
    (for xk in roots)
    (with n = (1- (length roots)))
    (collect (/ (* 2 (- 1 (expt xk 2)))
                (* (expt (+ n 1) 2)
                   (expt (funcall (legendre-0 n) xk) 2))))))

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

(defmacro define-integrator (name fn var int-value)
  (alexandria:with-gensyms (gint gfn)
    `(defmacro ,name (supports)
       (let ((,gint ,int-value))
         (flet ((,gfn (,var) ,fn))
           `(list ,@(mapcar (lambda (n)
                              (let* ((int (funcall (gauss-quadratur n) #',gfn))
                                     (err (abs (- int ,gint))))
                                `(cons ,int ,err)))
                            supports)))))))

(define-integrator herwig (log (+ x 2)) x (- (log 27.0d0) 2))
(define-integrator joe (/ 1 (+ 2 x)) x (log 3.0d0))

(defmacro main (int n-values)
  (alexandria:with-gensyms (gn gresult gint gerr)
    `(iter
       (with ,gresult = (time (,int ,n-values)))
       (for ,gn in (list ,@n-values))
       (for (,gint . ,gerr) in ,gresult)
       (format t "n: ~2,,d     int: ~,16E     err: ~,16E~%" ,gn ,gint ,gerr))))

(main herwig (2 4 8 16 32))
(main joe (2 4 8 16 32))
