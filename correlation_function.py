import numpy

def yaracorrelate(a, b, method='full'):
    a = (a - numpy.mean(a)) / (numpy.std(a) * len(a))
    b = (b - numpy.mean(b)) / (numpy.std(b))
    return numpy.correlate(a, b, method)


def yaracorrelate2(a, b, method='full', sps=1, onlypos=False):
    a = (a - numpy.mean(a)) / (numpy.std(a) * len(a))
    b = (b - numpy.mean(b)) / (numpy.std(b))
    corr = numpy.correlate(a, b, method)
    if onlypos:
        index_max = numpy.argmax(corr)
    else:
        index_max = numpy.argmax(abs(corr))
    if method == 'full':
        lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), (len(b) + len(a)) - 1) / sps
    elif method == 'same':
        lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), max(len(a), len(b))) / sps
    elif method == 'valid':
        lag = numpy.linspace(-(len(a) - 1), (len(b) - 1), (max(len(a), len(b)) - min(len(a), len(b)) + 1))

    lag_part = lag[index_max]
    corr_part = corr[index_max]
    return lag, corr, lag_part, corr_part


def corr_coeff(data1='', data2='', maxf='', minf=''):
    n = len(data1)
    n2 = len(data2)
    if n != n2:
        print('Datasets do not have same length!')
        return
    else:
        if max == '':
            r1 = n * numpy.sum(data1 * data2) - (numpy.sum(data1) * numpy.sum(data2))
            r2 = numpy.sqrt((n * numpy.sum(data1 ** 2) - numpy.sum(data1) ** 2) * (
                        n * numpy.sum(data2 ** 2) - numpy.sum(data2) ** 2))
            r = r1 / r2
            print('correlation coefficient r = ', r)
            return r
        if maxf != '':
            newdata2 = []
            newdata1 = []
            for i in range(0, n2):
                if data2[i] < maxf:
                    if data2[i] > minf:
                        newdata2.append(data2[i])
                        newdata1.append(data1[i])
                    else:
                        continue
                else:
                    continue
            data1 = numpy.asarray(newdata1)
            data2 = numpy.asarray(newdata2)
            r1 = n * numpy.sum(data1 * data2) - (numpy.sum(data1) * numpy.sum(data2))
            r2 = numpy.sqrt((n * numpy.sum(data1 ** 2) - numpy.sum(data1) ** 2) * (
                        n * numpy.sum(data2 ** 2) - numpy.sum(data2) ** 2))
            r = r1 / r2
            print('correlation coefficient r = ', r)
            return r