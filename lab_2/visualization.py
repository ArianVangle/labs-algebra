import matplotlib.pyplot as plt
import numpy as np


def plot_misclassified(X, y, y_pred, save_path="plots/misclassified.png"):
    errors = (y != y_pred)
    plt.figure(figsize=(6,5))
    plt.scatter(X[~errors,0], X[~errors,1], c=y[~errors], cmap='viridis', label='Correct')
    plt.scatter(X[errors,0], X[errors,1], c='red', marker='x', s=100, label='Misclassified')
    plt.title("Ошибки классификации"); plt.legend(); plt.xlabel("x1"); plt.ylabel("x2")
    plt.savefig(save_path, dpi=150); plt.close()


def plot_learning_curve(train_l, val_l, title, save):
    plt.figure(figsize=(6,4))
    plt.plot(train_l, label='Train Loss')
    plt.plot(val_l, label='Val Loss')
    plt.title(title); plt.xlabel("Epoch"); plt.ylabel("Loss"); plt.legend(); plt.grid()
    plt.savefig(save, dpi=150); plt.close()

def plot_decision_boundary(model, X, y, title, save):
    x_min, x_max = X[:,0].min()-0.5, X[:,0].max()+0.5
    y_min, y_max = X[:,1].min()-0.5, X[:,1].max()+0.5
    xx, yy = np.meshgrid(np.linspace(x_min, x_max, 200), np.linspace(y_min, y_max, 200))
    Z = model.predict_proba(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape)
    plt.figure(figsize=(6,5))
    plt.contourf(xx, yy, Z, alpha=0.3, cmap='RdBu')
    plt.scatter(X[:,0], X[:,1], c=y, cmap='RdBu', edgecolor='k', s=30)
    plt.title(title); plt.xlabel("x1"); plt.ylabel("x2"); plt.colorbar()
    plt.savefig(save, dpi=150); plt.close()
