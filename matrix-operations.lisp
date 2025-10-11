;;; Matrix Utilities in Common Lisp
;;; --------------------------------
;;; This file provides basic matrix operations (creation, addition,
;;; multiplication, transpose, row operations), and Gauss–Jordan elimination
;;; for solving systems of linear equations.

;;; --------------------------
;;; Matrix construction basics
;;; --------------------------

;; Create a new matrix of size rows x cols, with all entries initialized
;; to INIT (default 0).
(defun make-matrix (rows cols &optional (init 0))
  (make-array (list rows cols) :initial-element init))

;; Print a matrix in a nice row/column format.
(defun print-matrix (M)
  (dotimes (i (array-dimension M 0)) ; iterate over rows
    (dotimes (j (array-dimension M 1)) ; iterate over columns
      (format t "~3a " (aref M i j)))  ; print each element with width 3
    (terpri)))                         ; newline at the end of each row

;; Build an identity matrix of size n x n.
(defun identity-matrix (n)
  (let ((m (make-matrix n n)))
    (dotimes (i n)
      (setf (aref m i i) 1))  ; set diagonal entries to 1
    m))

;;; --------------------------
;;; Matrix arithmetic
;;; --------------------------

;; Add two matrices elementwise. Dimensions must match.
(defun matrix-add (A B)
  (unless (and (= (array-dimension A 0) (array-dimension B 0))
               (= (array-dimension A 1) (array-dimension B 1)))
    (error "Matrix dimensions do not match: ~A vs ~A"
           (array-dimensions A) (array-dimensions B)))
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (C (make-matrix rows cols)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref C i j) (+ (aref A i j) (aref B i j)))))
    C))

;; Multiply two matrices A and B. Number of columns in A must match
;; number of rows in B.
(defun matrix-multiply (A B)
  (unless (= (array-dimension A 1) (array-dimension B 0))
    (error "Matrix dimensions do not allow multiplication: ~A vs ~A"
           (array-dimensions A) (array-dimensions B)))
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension B 1))
         (inner (array-dimension A 1))
         (C (make-matrix rows cols 0)))
    (dotimes (i rows)
      (dotimes (j cols)
        (dotimes (k inner)
          (incf (aref C i j) (* (aref A i k) (aref B k j))))))
    C))

;; Return the transpose of matrix M (swap rows/columns).
(defun transpose (M)
  (let* ((rows (array-dimension M 0))
         (cols (array-dimension M 1))
         (MT (make-matrix cols rows)))
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref MT j i) (aref M i j))))
    MT))

;;; --------------------------
;;; Row operations
;;; --------------------------

;; Swap two rows (i and j) in place. Rotatef avoids a temp variable.
(defun swap-rows (M i j)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (rotatef (aref M i k) (aref M j k)))))

;; Multiply an entire row by a scalar factor in place.
(defun scale-row (M i factor)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (setf (aref M i k) (* factor (aref M i k))))))

;; Add factor * (row src) to (row dest). src is unchanged.
(defun add-rows (M src dest factor)
  (let ((cols (array-dimension M 1)))
    (dotimes (k cols)
      (incf (aref M dest k) (* factor (aref M src k))))))

;;; --------------------------
;;; Gauss–Jordan elimination
;;; --------------------------

;; Reduce augmented matrix M to Reduced Row Echelon Form (RREF).
;; Works in place, returns M.
(defun gauss-jordan (M)
  (let* ((rows (array-dimension M 0))
         (cols (array-dimension M 1))
         (lead 0))
    (dotimes (r rows M)  ; loop over rows
      (when (>= lead cols)
        (return M))      ; stop if we run out of columns

      ;; Find a pivot row
      (let ((i r))
        (loop while (and (< i rows) (zerop (aref M i lead)))
              do (incf i))
        (if (= i rows)
            ;; no pivot in this column, move right
            (progn (incf lead) (decf r)) 
            (progn
              ;; swap pivot row into place
              (swap-rows M r i)
              ;; normalize pivot row so pivot = 1
              (let ((pivot (aref M r lead)))
                (when (zerop pivot)
                  (error "Singular matrix, cannot invert/solve"))
                (scale-row M r (/ 1.0 pivot)))
              ;; eliminate pivot from all other rows
              (dotimes (j rows)
                (unless (= j r)
                  (let ((factor (aref M j lead)))
                    (when (/= factor 0)
                      (add-rows M r j (- factor))))))
              (incf lead)))))))

;; Solve system A x = b using Gauss–Jordan elimination.
;; Builds augmented matrix [A|b], reduces it, then extracts x.
(defun solve-system (A b)
  (let* ((rows (array-dimension A 0))
         (cols (array-dimension A 1))
         (M (make-matrix rows (1+ cols))))
    ;; build augmented matrix
    (dotimes (i rows)
      (dotimes (j cols)
        (setf (aref M i j) (aref A i j)))
      (setf (aref M i cols) (aref b i 0)))
    ;; row reduce
    (gauss-jordan M)
    ;; extract solution vector x
    (let ((x (make-matrix cols 1)))
      (dotimes (i cols)
        (setf (aref x i 0) (aref M i cols)))
      x)))
