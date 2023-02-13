import numpy as np
from scipy import sparse
from scipy.sparse.linalg import lsmr
def GenInv(mode, logO, logR, iter, *arg):
    '''
    Generalized Inversion (Andrews 1986)
    equation for each frequency f: Oij = Ei * Sj * (Rij)^-1 * exp(-(pi*Rij*f)/(Qs*Vs))
    equation after taking log10: log10Oij + log10Rij = log10Ei + log10Sj - (pi*Rij*f*log10(e))/Vs * Qs^-1
    model parameters: log10Ei, log10Sj, Qs^-1
    There are 6 ways of inversion:
        mode = 1: path (geometric), only source term, no site term
        mode = 2: path (geometric), both source and Q term, no site term
        mode = 3: path (geometric), both source and site term
        mode = 4: path (geometric), all source, site and Q term
        mode = 5: path (geometric), then only source term, site term from input (e.g. synthetic)
        mode = 6: path (geometric), then both source term and Q, site term from input (e.g. synthetic)
        ** note that for mode 1 and 5, no need for inversion, just take mean from all channel measurements

    input: (all of the terms are log10(amplitudes))
    :param mode: 1 to 6
    :param logO: data Oij (size of nev * nch * nf)
    :param logR: path term Rij (size of nev * nch)
    :param iter: maximum iterative time of reweighted lsqr
    :param arg:
        if mode == 1, none
        if mode == 2, constant Vs
        if mode == 3, none
        if mode == 4, constant Vs
        if mode == 5, site term logS (size of nch * nf)
        if mode == 6, constant Vs, site term logS (size of nch * nf)
    :return:
    logE: (size of nev * nf)
    logS: (size of nch * nf)
    Qinv: (size of nf)
    syndata: same shape as logO
    '''



    # t1 = time()
    eps = 1e-1
    nf = logO.shape[2]
    nev = logO.shape[0]
    nch = logO.shape[1]
    logE = np.zeros((nev, nf))
    logE[:] = np.nan
    logS = np.zeros((nch, nf))
    logS[:] = np.nan
    Qinv = np.zeros(nf)
    Qinv[:] = np.nan
    syndata = np.zeros_like(logO)
    syndata[:] = np.nan
    weights = np.ones_like(logO)

    if mode == 1:
        for ii in range(iter):
            for idx in range(nf):
                logE[:, idx] = np.sum((logO[:, :, idx] + logR) * weights[:,:,idx], 1) / np.sum(weights[:,:,idx], 1)

            syndata = np.repeat(np.expand_dims(logE, axis=1), nch, 1) - np.repeat(np.expand_dims(logR, axis=2), nf, 2)
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)

    elif mode == 2:
        vs = arg[0]
        row = np.tile(np.arange(nch * nev), 2)
        col = np.concatenate((np.repeat(np.arange(nev), nch), nev * np.ones(nev * nch)))
        val = np.ones(nch * nev )
        val = np.concatenate((val, -np.pi * np.log10(np.e) / vs * (10**logR).reshape(-1)))
        G = sparse.coo_matrix((val, (row, col)), shape=(nch * nev, nev + 1))
        for ii in range(iter):
            for idx in range(nf):
                if np.isnan(logO[0, 2, idx]) < 1:
                    W = sparse.diags(weights[:, :, idx].reshape(-1))
                    m = lsmr(W.dot(G), W.dot((logO[:, :, idx] + logR).reshape(-1)))[0]
                    logE[:, idx] = m[:nev] * 1
                    Qinv[idx] = m[-1] * 1
                    syndata[:, :, idx] = G.dot(m).reshape(nev, nch) - logR
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)

    elif mode == 3:
        row = np.tile(np.arange(nch * nev), 2)
        col = np.concatenate((np.repeat(np.arange(nev), nch), np.tile(nev + np.arange(nch), nev)))
        val = np.ones(nch * nev * 2)
        G = sparse.coo_matrix((val, (row, col)), shape=(nch * nev, nch + nev))
        for ii in range(iter):
            for idx in range(nf):
                if np.isnan(logO[0, 2, idx]) < 1:
                    W = sparse.diags(weights[:, :, idx].reshape(-1))
                    m = lsmr(W.dot(G), W.dot((logO[:, :, idx] + logR).reshape(-1)))[0]
                    logE[:, idx] = m[:nev] *1
                    logS[:, idx] = m[nev:]*1
                    syndata[:, :, idx] = G.dot(m).reshape(nev, nch) - logR
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)

    elif mode == 4:
        vs = arg[0]
        row = np.tile(np.arange(nch * nev), 2)
        col = np.concatenate((np.repeat(np.arange(nev), nch), np.tile(nev + np.arange(nch), nev)))
        val = np.ones(nch * nev * 2)
        row = np.concatenate((row, np.arange(nev * nch)))
        col = np.concatenate((col, (nch + nev) * np.ones(nev * nch)))
        val = np.concatenate((val, -np.pi * np.log10(np.e) / vs * (10 ** logR).reshape(-1)))
        G = sparse.coo_matrix((val, (row, col)), shape=(nch * nev, nch + nev + 1))
        for ii in range(iter):
            for idx in range(nf):
                if np.isnan(logO[0, 2, idx]) < 1:
                    W = sparse.diags(weights[:, :, idx].reshape(-1))
                    m = lsmr(W.dot(G), W.dot((logO[:, :, idx] + logR).reshape(-1)))[0]
                    logE[:, idx] = m[:nev]*1
                    logS[:, idx] = m[nev:nev + nch]*1
                    Qinv[idx] = m[-1]*1
                    syndata[:, :, idx] = G.dot(m).reshape(nev, nch) - logR
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)


    elif mode == 5:
        logS = arg[0]
        for idx in range(nf):
            logE[:, idx] = np.nanmean(logO[:,:,idx]+logR-np.repeat(logS[:,idx].reshape(1,nch),nev,0), 1)
        syndata = np.repeat(np.expand_dims(logE, axis=1), nch, 1
                            ) - np.repeat(np.expand_dims(logR, axis=2), nf, 2
                                          ) + np.repeat(np.expand_dims(logS, axis=0),nev,0)
        for ii in range(iter):
            for idx in range(nf):
                logE[:, idx] = np.sum((logO[:, :, idx] + logR-np.repeat(logS[:,idx].reshape(1,nch),nev,0)
                                       ) * weights[:,:,idx], 1) / np.sum(weights[:,:,idx], 1)
            syndata = np.repeat(np.expand_dims(logE, axis=1), nch, 1
                                    ) - np.repeat(np.expand_dims(logR, axis=2), nf, 2
                                                  ) + np.repeat(np.expand_dims(logS, axis=0), nev, 0)
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)


    elif mode == 6:
        vs = arg[0]
        logS = arg[1]
        row = np.tile(np.arange(nch * nev), 2)
        col = np.concatenate((np.repeat(np.arange(nev), nch), nev * np.ones(nev * nch)))
        val = np.ones(nch * nev)
        val = np.concatenate((val, -np.pi * np.log10(np.e) / vs * (10 ** logR).reshape(-1)))
        G = sparse.coo_matrix((val, (row, col)), shape=(nch * nev, nev + 1))

        for ii in range(iter):
            for idx in range(nf):
                if np.isnan(logO[0, 2, idx]) < 1:
                    W = sparse.diags(weights[:, :, idx].reshape(-1))
                    m = lsmr(W.dot(G), W.dot((logO[:, :, idx] + logR-np.repeat(logS[:,idx].reshape(1,nch),nev,0
                                                                               )).reshape(-1)))[0]
                    logE[:, idx] = m[:nev] * 1
                    Qinv[idx] = m[-1] * 1
                    syndata[:, :, idx] = G.dot(m).reshape(nev, nch) - logR+ np.repeat(logS[:,idx].reshape(1,nch),nev,0)
            weights = 1 / (eps + np.abs(logO - syndata) ** 1)

    else:
        print('''GenInv(mode, frequency, logO, logR, *arg), 
        There are 6 ways of inversion:
        mode = 1: path (geometric), only source term, no site term
        mode = 2: path (geometric), both source and Q term, no site term
        mode = 3: path (geometric), both source and site term
        mode = 4: path (geometric), all source, site and Q term
        mode = 5: path (geometric), then only source term, site term from input (e.g. synthetic)
        mode = 6: path (geometric), then both source term and Q, site term from input (e.g. synthetic)
        \n''')

    # print(time() - t1)
    return logE, logS, Qinv, syndata
