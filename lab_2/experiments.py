import numpy as np
import matplotlib.pyplot as plt
from lab_2.model import SingleLayerPerceptron

def run_param_experiment(name, configs, X_tr, X_te, y_tr, y_te, X_v, y_v, param_type="lr"):
    print(f"🔹 Эксперимент: {name}")
    plt.figure(figsize=(6,4))
    for val, label in configs:
        params = {"w_init": "random_small"}
        fit_params = {"epochs": 100, "batch_size": 32, "lr": 0.1}
        
        if param_type == "w_init": params["w_init"] = val
        elif param_type == "batch_size": fit_params["batch_size"] = val
        elif param_type == "lr": fit_params["lr"] = val

        m = SingleLayerPerceptron(**params)
        m.fit(X_tr, y_tr, X_v, y_v, **fit_params)
        plt.plot(m.train_loss_hist, label=label)
        
    plt.title(f"Влияние {param_type}"); plt.legend(); plt.grid()
    plt.savefig(f"plots/exp_{name.replace(' ','_')}.png"); plt.close()

def cross_validate(X, y, lr_grid, bs_grid, k=5):
    fold_size = len(y) // k
    best_score, best_cfg = 0, None
    
    for lr in lr_grid:
        for bs in bs_grid:
            scores = []
            for i in range(k):
                val_idx = np.arange(i*fold_size, (i+1)*fold_size)
                tr_idx = np.delete(np.arange(len(y)), val_idx)
                
                m = SingleLayerPerceptron()
                m.fit(X[tr_idx], y[tr_idx], X[val_idx], y[val_idx], epochs=50, lr=lr, batch_size=bs)
                scores.append(np.mean(m.predict(X[val_idx]) == y[val_idx]))
            
            avg = np.mean(scores)
            if avg > best_score:
                best_score, best_cfg = avg, (lr, bs)
    return best_cfg, best_score


def run_experiment(name, configs, X_tr, X_te, y_tr, y_te, X_v, y_v, param_name="lr"):
    """Универсальный раннер экспериментов"""
    print(f"\n🔹 Эксперимент: {name}")
    results = []
    plt.figure(figsize=(6,4))
    for val, label in configs:
        init_kwargs, fit_kwargs = {}, {"epochs": 100, "batch_size": 32, "lr": 0.1}
        if param_name == "w_init": init_kwargs["w_init"] = val
        elif param_name == "batch_size": fit_kwargs["batch_size"] = val
        elif param_name == "lr": fit_kwargs["lr"] = val

        m = SingleLayerPerceptron(**init_kwargs)
        m.fit(X_tr, y_tr, X_v, y_v, **fit_kwargs)
        acc_tr = np.mean(m.predict(X_tr) == y_tr)
        acc_te = np.mean(m.predict(X_te) == y_te)
        results.append(f"{label}: Train={acc_tr:.3f}, Test={acc_te:.3f}")
        plt.plot(m.train_loss_hist, label=label)
    plt.title(f"Влияние {param_name}"); plt.xlabel("Epoch"); plt.ylabel("Loss"); plt.legend(); plt.grid()
    plt.savefig(f"plots/exp_{name.replace(' ','_')}.png", dpi=150); plt.close()
    print("\n".join(results))

def run_kfold_cv(X, y, lr_grid, bs_grid, k=5):
    print("\n🔹 Задание 5: 5-Fold Cross-Validation")
    fold_size = len(y) // k
    best_score, best_cfg = 0, {}
    for lr in lr_grid:
        for bs in bs_grid:
            scores = []
            for i in range(k):
                val_idx = np.arange(i*fold_size, (i+1)*fold_size)
                tr_idx = np.concatenate([np.arange(0, i*fold_size), np.arange((i+1)*fold_size, len(y))])
                X_tr, X_val = X[tr_idx], X[val_idx]
                y_tr, y_val = y[tr_idx], y[val_idx]
                m = SingleLayerPerceptron()
                m.fit(X_tr, y_tr, X_val, y_val, epochs=50, lr=lr, batch_size=bs)
                scores.append(np.mean(m.predict(X_val) == y_val))
            avg, std = np.mean(scores), np.std(scores)
            if avg > best_score: best_score, best_cfg = avg, (lr, bs)
            print(f"  lr={lr}, bs={bs} -> Acc={avg:.4f} ± {std:.4f}")
    print(f"Лучшие параметры: lr={best_cfg[0]}, batch_size={best_cfg[1]} (Acc={best_score:.4f})")
    return best_cfg