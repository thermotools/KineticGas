def equal(a, b, tol=1e-10):
    return abs((a - b) / max(abs(a), abs(b))) < tol