import numpy as np

def standardize(X_train, X_test):
    mean = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1.0
    return (X_train - mean) / std, (X_test - mean) / std

def generate_custom_data(mode="gaussian", noise_prob=0.0, n=300, seed=42):
    rng = np.random.RandomState(seed)
    if mode == "gaussian":
        X1 = rng.multivariate_normal([2, 2], [[1, 0.5], [0.5, 1]], n//2)
        X2 = rng.multivariate_normal([-2, -2], [[1, -0.5], [-0.5, 1]], n//2)
        X = np.vstack([X1, X2])
        y = np.hstack([np.zeros(n//2), np.ones(n//2)])
    elif mode == "circle":
        R = rng.rand(n) * 2 + 0.5
        Theta = rng.rand(n) * 2 * np.pi
        X = np.c_[R * np.cos(Theta), R * np.sin(Theta)]
        y = (R > 2.0).astype(int)
    elif mode == "xor":
        X = np.array([[0,0], [0,1], [1,0], [1,1]] * (n//4))
        y = np.array([0, 1, 1, 0] * (n//4))
        
    if noise_prob > 0:
        flip_mask = rng.rand(len(y)) < noise_prob
        y[flip_mask] = 1 - y[flip_mask]
    return X, y

def compute_metrics(y_true, y_pred, y_proba):
    tp = np.sum((y_true == 1) & (y_pred == 1))
    tn = np.sum((y_true == 0) & (y_pred == 0))
    fp = np.sum((y_true == 0) & (y_pred == 1))
    fn = np.sum((y_true == 1) & (y_pred == 0))
    acc = (tp+tn)/(tp+tn+fp+fn) if (tp+tn+fp+fn)>0 else 0
    prec = tp/(tp+fp) if (tp+fp)>0 else 0
    rec = tp/(tp+fn) if (tp+fn)>0 else 0
    f1 = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0
    
    thresholds = np.sort(np.unique(y_proba))[::-1]
    tpr_list, fpr_list = [], []
    pos = np.sum(y_true == 1)
    neg = np.sum(y_true == 0)
    for thr in thresholds:
        p = (y_proba >= thr).astype(int)
        tpr_list.append(np.sum((y_true==1)&(p==1))/pos if pos>0 else 0)
        fpr_list.append(np.sum((y_true==0)&(p==1))/neg if neg>0 else 0)
    tpr_arr, fpr_arr = np.array(tpr_list), np.array(fpr_list)
    auc = np.sum((fpr_arr[1:] - fpr_arr[:-1]) * (tpr_arr[1:] + tpr_arr[:-1]) / 2)
    return acc, prec, rec, f1, auc, tpr_arr, fpr_arr