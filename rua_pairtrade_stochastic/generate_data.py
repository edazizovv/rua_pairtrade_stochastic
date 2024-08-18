import numpy


def handle_generator(M_x0, M_r, B_x0, B_mu, B_sigma, X_x0, X_k, X_theta, X_eta, n):
    # generate dM(t) and M(t)
    dM, M = riskfreeasset(x0=M_x0, r=M_r, n=n)

    # generate dB(t) and B(t)
    dB, B = geometricBrownian(x0=B_x0, mu=B_mu, sigma=B_sigma, n=n)

    # generate dX(t) and X(t)
    dX, X = generateOrsteinUhlenbeck(x0=X_x0, k=X_k, theta=X_theta, eta=X_eta, n=n)

    # calculate dA(t) and A(t)
    A = numpy.exp(X) * B
    dA = numpy.diff(A)

    # return dM(t), M(t), dA(t), A(t), dB(t), B(t), dX(t), X(t)
    return dM, M, dA, A, dB, B, dX, X


def riskfreeasset(x0, r, n):
    # sampling
    process_increments = numpy.array([r] * n, dtype=numpy.float64)
    cum_prod = numpy.cumprod(a=process_increments, axis=0, dtype=numpy.float64)
    process_values = numpy.concatenate((numpy.array([x0], dtype=numpy.float64), x0 * cum_prod), axis=0)

    # return dM(t), M(t)
    return process_increments, process_values


def geometricBrownian(x0, mu, sigma, n):
    # generate normals
    standard_normal_sample = numpy.random.normal(loc=0.0, scale=1.0, size=(n,))

    # compute process
    process_increments = numpy.array([], dtype=numpy.float64)
    process_values = numpy.array([x0], dtype=numpy.float64)
    for t in numpy.arange(n):
        # an increment
        an_increment = mu * process_values[t] + sigma * standard_normal_sample[t]
        # a value
        a_value = process_values[t] + an_increment
        # arrays' update
        process_increments = numpy.append(arr=process_increments, values=numpy.array([an_increment], dtype=numpy.float64), axis=0)
        process_values = numpy.append(arr=process_values, values=numpy.array([a_value], dtype=numpy.float64), axis=0)

    # return dZ(t), Z(t)
    return process_increments, process_values


def generateOrsteinUhlenbeck(x0, k, theta, eta, n):
    # generate normals
    standard_normal_sample = numpy.random.normal(loc=0.0, scale=1.0, size=(n,))

    # compute process
    process_increments = numpy.array([], dtype=numpy.float64)
    process_values = numpy.array([x0], dtype=numpy.float64)
    for t in numpy.arange(n):
        # an increment
        an_increment = k * (theta - process_values[t]) + eta * standard_normal_sample[t]
        # a value
        a_value = process_values[t] + an_increment
        # arrays' update
        process_increments = numpy.append(arr=process_increments, values=numpy.array([an_increment], dtype=numpy.float64), axis=0)
        process_values = numpy.append(arr=process_values, values=numpy.array([a_value], dtype=numpy.float64), axis=0)

    # return dZ(t), Z(t)
    return process_increments, process_values
