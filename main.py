#----------------------
# Shay Fletcher     318727641
# Nika Tatikishvili 321735433
# https://github.com/shayfletcherz/NumericalAnalysisProject5.git
#----------------------


def linear_calc(x1, y1, x2, y2, xf):
    return ((y1 - y2) / (x1 - x2)) * xf + ((y2 * x1) - (y1 * x2)) / (x1 - x2)



def linear_interpolation(table, xf):
    for row in range(len(table) - 1):
        if (xf > table[row][0]) and xf < table[row + 1][0]:
            x1 = table[row][0]
            x2 = table[row + 1][0]
            y1 = table[row][1]
            y2 = table[row + 1][1]
            return (((y1 - y2) / (x1 - x2)) * xf) + ((y2 * x1 - y1 * x2) / (x1 - x2))


def lagrange_interpolation(table, xf):
    y = 0
    for k in range(len(table)):
        if len(table) != len(table[k]):
            print("ERROR")
            return 1
        t = 1
        for j in range(len(table)):
            if j != k:
                t = t * ((xf - table[j]) / (table[k] - table[j]))
        y += t * table[k][1]
    return y

def neville_interpolation(table,xf):
    n = len(table)
    x = 0
    y = 1
    for i in range(1,n,+1):
        for j in range(n-1,i-1,-1):
            table[j][y] = ((xf-table[j-i][x])*table[j][y]-(xf-table[j][x])*table[j-1][y])/(table[j][x]-table[j-i][x])
    result = table[n-1][y]
    return result

#Creating a polynomial matrix
def makePolynomialMat(points):
    size = len(points)
    newMat = makeMatrics(size, size)
    newB = makeMatrics(size, size)
    for i in range(size):
        xi = points[i][0]
        for j in range(size):
            newMat[i][j] = xi ** j
        newB[i][0] = points[i][1]
    return newMat, newB


#Copy matrix function
def copyMat(matrix):
    B = makeMatrics(len(matrix), len(matrix[0]))
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            B[i][j] = matrix[i][j]
    return B

#Exchanging cols
def manualSwapCol(matrix, col1, col2):
    if col2 < len(matrix) and col1 < len(matrix):
        for i in range(len(matrix)):
            temp = matrix[i][col1]
            matrix[i][col1] = matrix[i][col2]
            matrix[i][col2] = temp
    return matrix


#Exchanging lines
def manualSwapRow(matrix, b, row1, row2):
    if row2 < len(matrix) and row1 < len(matrix):
        temp = matrix[row1]
        matrix[row1] = matrix[row2]
        matrix[row2] = temp
        if b is not None:
            temp = b[row1]
            b[row1] = b[row2]
            b[row2] = temp
    return matrix, b


def rowSum(line):
    lineSum = 0
    for index in range(len(line)):  # run over all the line`s members
        lineSum += abs(line[index])
    return lineSum


def createDominantDiagonal(matrix, b):
    max = 0
    maxIndex = 0
    sum = 0
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            sum += abs(matrix[i][j])
            if abs(matrix[i][j]) > max:
                max = abs(matrix[i][j])
                maxIndex = j
        if (sum - max) <= max:
            matrix = manualSwapCol(matrix, maxIndex, i)
        else:
            max = 0
            maxIndex = 0
            for j in range(len(matrix)):
                sum += abs(matrix[j][i])
                if abs(matrix[j][i]) > max:
                    max = abs(matrix[j][i])
                    maxIndex = j
            if rowSum(matrix[j]) - max <= max:
                matrix, b = manualSwapRow(matrix, b, i, maxIndex)
            else:
                return None, None
    return matrix, b


#Creating a zero matrix
def makeMatrics(row, col=1):
    c = []
    for i in range(row):
        c += [[0] * col]
    return c


#Presentation matrix and vector result after pivot
def isolateVariables(coefficientMat, b):
    vectorB = makeMatrics(len(coefficientMat))
    matA = makeMatrics(len(coefficientMat), len(coefficientMat[0]))
    for i in range(len(coefficientMat)):
        vectorB[i][0] = b[i][0] / coefficientMat[i][i]
        j = 0
        while j < len(coefficientMat[0]):
            if j is i:
                matA[i][i] = 1
            else:
                matA[i][j] -= coefficientMat[i][j] / coefficientMat[i][i]
            j += 1
    return matA, vectorB


#Presentation of iterations using the Seidel method
def gaussSeidelIter(coefficientMat, b):
    epsilon = 0.00001
    iteration = 0
    flag = True
    prevX = makeMatrics(len(coefficientMat))  # start as zero vector
    currentX = makeMatrics(len(coefficientMat))
    matA, vectorB = isolateVariables(coefficientMat, b)
    while abs(currentX[0][0] - prevX[0][0]) > epsilon or flag is True:
        flag = False
        prevX[0][0] = currentX[0][0]
        if iteration >= 100:
            break
        for i in range(len(coefficientMat)):
            j = 0
            currentX[i][0] = vectorB[i][0]
            while j < len(coefficientMat[0]):
                if j is not i:
                    currentX[i][0] += matA[i][j] * currentX[j][0]
                j += 1
        iteration += 1
    return currentX


# A function will get coefficients of the polynomial and the X that we want to find and return an approximate value
def getCoefficientsCalcY(coefficients, X):
    sum = 0
    for x in range(len(coefficients)):
        sum += coefficients[x][0] * (X ** x)
    return sum


#Calculation of an estimated value by a polynomial method
def polynomial(table, X):
    a, b = makePolynomialMat(table)
    copyA = copyMat(a)
    copyB = copyMat(b)
    copyA, copyB = createDominantDiagonal(copyA, copyB)
    if (copyA is not None) and (copyB is not None):
        a = copyA
        b = copyB
    coefficients = gaussSeidelIter(a, b)
    print("The Formula: y(f)= y1-y2        y2*x1-y1*x2")
    print("                  -----  * xf + -----------")
    print("                   x1-x2         x1-x2\n")
    print("The Approximate value by Polynomial interpolation:")
    x = str(getCoefficientsCalcY(coefficients, X))
    return x


def main():
    user_choice = int((input("Please choose an option:\n 1 For Neville Interpolation,\n 2 For Polynomial Interpolation,\n 3 For Linear Interpolation,\n 4 For Lagrange Interpolation,\n 5 For Everything And Anything Else is Exit \n")))
    table = [[1.2, 3.5095],
             [1.3, 3.6984],
             [1.4, 3.9043],
             [1.5, 4.1293],
             [1.6, 4.3756]]
    point = 1.37
    linPoint = 3
    lagPoint = 1.43
    print('The points are:', table)
    print('The point:', point)

    if user_choice == 1:
        print("Neville interpolation:")
        print(neville_interpolation(table, point))

    elif user_choice == 2:
        print("Polynomial interpolation:")
        print(polynomial(table, point))

    elif user_choice == 3:
        print("Linear interpolation:")
        print(linear_interpolation(table, point))

    elif user_choice == 4:
        print("Lagrange interpolation:")
        print(lagrange_interpolation(table, lagPoint))

    elif user_choice == 5:
        print("Neville interpolation:")
        print(neville_interpolation(table, point))
        print("\n")
        print("Polynomial interpolation:")
        print(polynomial(table, point))
        print("\n")
        print("Linear interpolation:")
        print(linear_interpolation(table, linPoint))
        print("\n")
        print("Lagrange interpolation:")
        print(lagrange_interpolation(table, point))

    else:
        print("GoodBye")

main()