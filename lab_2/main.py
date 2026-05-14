import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
import json
import os

from lab_2.model import SingleLayerPerceptron
from lab_2.utils import *
from lab_2.visualization import *
from lab_2.experiments import *

os.makedirs("plots", exist_ok=True)



if __name__ == "__main__":
    print("="*60 + "\nЗАПУСК ЛАБОРАТОРНОЙ РАБОТЫ\n" + "="*60)

    print("\nГенерация данных (make_classification, n=500, 2 признака)...")
    X, y = make_classification(n_samples=500, n_features=2, n_redundant=0, n_informative=2, 
                               n_clusters_per_class=1, random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, stratify=y_train, random_state=42)
    
    X_train_s, X_test_s = standardize(X_train, X_test)
    std_tr = X_train.std(axis=0); std_tr[std_tr==0]=1.0
    X_val_s = (X_val - X_train.mean(axis=0)) / std_tr

    print("\nОбучение...")
    model = SingleLayerPerceptron()
    model.fit(X_train_s, y_train, X_val_s, y_val, epochs=100, lr=0.1, batch_size=32)
    
    plot_learning_curve(model.train_loss_hist, model.val_loss_hist, "Loss (Train vs Val)", "plots/loss_curve.png")
    plot_decision_boundary(model, X_test_s, y_test, "Decision Boundary", "plots/boundary.png")
    
    acc_tr = np.mean(model.predict(X_train_s) == y_train)
    acc_te = np.mean(model.predict(X_test_s) == y_test)
    print(f"Accuracy: Train={acc_tr:.4f}, Test={acc_te:.4f}")

    run_experiment("Learning Rate", 
                   [(0.001,"lr=0.001"), (0.01,"lr=0.01"), (0.5,"lr=0.5"), (1.0,"lr=1.0")],
                   X_train_s, X_test_s, y_train, y_test, X_val_s, y_val, param_name="lr")
    run_experiment("Batch Size",
                   [(1,"bs=1"), (16,"bs=16"), (64,"bs=64"), (256,"bs=256")],
                   X_train_s, X_test_s, y_train, y_test, X_val_s, y_val, param_name="batch_size")
    run_experiment("Weight Init",
                   [("zero","w=0"), ("random_small","w~N(0,0.01)"), ("random_large","w~N(0,10)")],
                   X_train_s, X_test_s, y_train, y_test, X_val_s, y_val, param_name="w_init")

    print("\nДоп. 1: Кастомные данные и визуализация...")
    for mode in ["gaussian", "circle", "xor"]:
        Xc, yc = generate_custom_data(mode, noise_prob=0.05, n=400, seed=42)
        Xc_tr, Xc_te, yc_tr, yc_te = train_test_split(Xc, yc, test_size=0.3, stratify=yc, random_state=42)
        Xc_tr_s, Xc_te_s = standardize(Xc_tr, Xc_te)
    
        mc = SingleLayerPerceptron()
        mc.fit(Xc_tr_s, yc_tr, Xc_tr_s, yc_tr, epochs=100, lr=0.2, batch_size=32)
    
        acc = np.mean(mc.predict(Xc_te_s) == yc_te)
        conclusion = "Линейно разделимо (перцептрон справляется)" if mode == "gaussian" else "Нелинейно (однослойный перцептрон не проведёт идеальную границу)"
        print(f"  {mode.upper()}: Accuracy={acc:.3f} {conclusion}")
    
        plot_decision_boundary(mc, Xc_te_s, yc_te, f"Decision Boundary: {mode.upper()}", f"plots/boundary_{mode}.png")
    
    print("\nВывод по Заданию 1:")
    print("   Однослойный перцептрон строит только линейную разделяющую гиперплоскость w^Tx+b=0.")
    print("   Он успешно разделяет гауссовы облака, но принципиально не может описать окружность или XOR без дополнительных признаков или скрытых слоёв.")

    print("\nДоп. 2: L2-регуляризация & Hinge Loss")
    m_l2 = SingleLayerPerceptron()
    m_l2.fit(X_train_s, y_train, X_val_s, y_val, epochs=100, lr=0.1, batch_size=32, l2_lambda=0.5)
    print(f"  L2 Lambda=0.5 -> Train Loss: {m_l2.train_loss_hist[-1]:.4f}")
    plot_learning_curve(m_l2.train_loss_hist, m_l2.val_loss_hist, "Loss with L2(λ=0.5)", "plots/loss_l2.png")
    
    m_hinge = SingleLayerPerceptron()
    m_hinge.fit(X_train_s, y_train, X_val_s, y_val, epochs=100, lr=0.1, batch_size=32, loss_type="hinge")
    print(f"  Hinge Loss -> Train Loss: {m_hinge.train_loss_hist[-1]:.4f}")

    print("\nДоп. 3: Метрики, ROC, Анализ ошибок...")
    preds = model.predict(X_test_s)
    acc, prec, rec, f1, auc, tpr, fpr = compute_metrics(y_test, preds, model.predict_proba(X_test_s))
    print(f"  Acc={acc:.3f} | Prec={prec:.3f} | Rec={rec:.3f} | F1={f1:.3f} | ROC-AUC={auc:.3f}")
    plt.figure(figsize=(5,5))
    plt.plot(fpr, tpr, label=f'ROC (AUC={auc:.3f})')
    plt.plot([0,1], [0,1], 'k--')
    plt.xlabel("FPR"); plt.ylabel("TPR"); plt.legend(); plt.grid()
    plt.savefig("plots/roc_curve.png", dpi=150); plt.close()
    plot_misclassified(X_test_s, y_test, preds)

    print("\nДоп. 4: Momentum SGD (β=0.5, 0.9, 0.99)")
    plt.figure(figsize=(6,4))
    for beta in [0.5, 0.9, 0.99]:
        mm = SingleLayerPerceptron()
        mm.fit(X_train_s, y_train, X_val_s, y_val, epochs=50, lr=0.1, batch_size=32, use_momentum=True, momentum_beta=beta)
        plt.plot(mm.train_loss_hist, label=f"β={beta}")
    plt.title("Momentum SGD Loss"); plt.xlabel("Epoch"); plt.ylabel("Loss"); plt.legend(); plt.grid()
    plt.savefig("plots/momentum_loss.png", dpi=150); plt.close()

    best_cfg = run_kfold_cv(X_train_s, y_train, lr_grid=[0.01, 0.1, 0.5], bs_grid=[16, 32, 64])
    print("\n🔹 Финальное обучение на всех обучающих данных...")

    X_train_full = np.vstack([X_train_s, X_val_s])
    y_train_full = np.hstack([y_train, y_val])

    final_model = SingleLayerPerceptron()
    final_model.fit(
        X_train_full, y_train_full,
        X_val=X_val_s, y_val=y_val,
        epochs=100,
        lr=best_cfg[0],
        batch_size=best_cfg[1]
    )

    acc_final = np.mean(final_model.predict(X_test_s) == y_test)
    preds_final = final_model.predict(X_test_s)
    proba_final = final_model.predict_proba(X_test_s)

    acc, prec, rec, f1, auc, tpr, fpr = compute_metrics(y_test, preds_final, proba_final)

    print(f"\nФИНАЛЬНЫЕ РЕЗУЛЬТАТЫ (на независимом тесте):")
    print(f"   Accuracy : {acc:.4f}")
    print(f"   Precision: {prec:.4f}")
    print(f"   Recall   : {rec:.4f}")
    print(f"   F1-score : {f1:.4f}")
    print(f"   ROC-AUC  : {auc:.4f}")
    print(f"   Hyperparams: lr={best_cfg[0]}, batch_size={best_cfg[1]}")

    report = {
        "best_lr": best_cfg[0],
        "best_batch_size": best_cfg[1],
        "final_test_accuracy": acc,
        "final_test_metrics": {
            "Accuracy": acc, "Precision": prec, "Recall": rec, 
            "F1": f1, "ROC_AUC": auc
        },
        "weights": final_model.w.tolist(),
        "bias": float(final_model.b)
    }
    with open("final_report.json", "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print("\nВсе графики сохранены в папке /plots/")
    print("="*60 + "\nЛАБОРАТОРНАЯ РАБОТА ЗАВЕРШЕНА\n" + "="*60)