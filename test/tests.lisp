;;; --------------------------
;;; Test suite for matrix utilities
;;; --------------------------

(defpackage :matrix-operations/test
  (:use :cl :matrix-operations))
(in-package :matrix-operations/test)


(defun test-matrix-utils ()
  (format t "~%=== Matrix Utility Tests ===~%")

  ;; 1️⃣ Create two small 2×2 matrices
  (let* ((A (make-matrix 2 2))
         (B (make-matrix 2 2)))
    (setf (aref A 0 0) 1 (aref A 0 1) 2
          (aref A 1 0) 3 (aref A 1 1) 4)
    (setf (aref B 0 0) 5 (aref B 0 1) 6
          (aref B 1 0) 7 (aref B 1 1) 8)

    (format t "~%A = ~%") (print-matrix A)
    (format t "~%B = ~%") (print-matrix B)

    ;; 2️⃣ Add them
    (format t "~%A + B = ~%")
    (print-matrix (matrix-add A B))

    ;; 3️⃣ Multiply them
    (format t "~%A × B = ~%")
    (print-matrix (matrix-multiply A B))

    ;; 4️⃣ Transpose of A
    (format t "~%Transpose(A) = ~%")
    (print-matrix (transpose A))

    ;; 5️⃣ Solve simple system: A·x = b
    ;; A = [[2, 1], [5, 7]], b = [11, 13]
    (let* ((A2 (make-matrix 2 2))
           (b (make-matrix 2 1)))
      (setf (aref A2 0 0) 2 (aref A2 0 1) 1
            (aref A2 1 0) 5 (aref A2 1 1) 7)
      (setf (aref b 0 0) 11 (aref b 1 0) 13)
      (format t "~%Solving A·x = b for~%A = ~%")
      (print-matrix A2)
      (format t "b = ~%")
      (print-matrix b)
      (let ((x (solve-system A2 b)))
        (format t "~%Solution x = ~%")
        (print-matrix x)))))

;;; Run tests automatically
(test-matrix-utils)
