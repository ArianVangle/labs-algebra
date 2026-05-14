import numpy as np

class SingleLayerPerceptron:
    def __init__(self, w_init="random_small", b_init=0.0, seed=42):
        self.w_init = w_init
        self.b_init = b_init
        self.seed = seed
        self.w = None
        self.b = None
        self.train_loss_hist = []
        self.val_loss_hist = []

    def _init_params(self, n_features):
        rng = np.random.RandomState(self.seed)
        if self.w_init == "zero":
            self.w = np.zeros(n_features)
        elif self.w_init == "random_small":
            self.w = rng.randn(n_features) * 0.01
        elif self.w_init == "random_large":
            self.w = rng.randn(n_features) * 10.0
        self.b = float(self.b_init)

    @staticmethod
    def sigmoid(z):
        z = np.clip(z, -500, 500)
        return 1.0 / (1.0 + np.exp(-z))

    def forward(self, X):
        """Прямой проход"""
        return self.sigmoid(X @ self.w + self.b)

    def compute_loss(self, y_true, y_pred_or_z, l2_lambda=0.0, use_hinge=False, weights=None):
        """Функция потерь"""
        eps = 1e-12
        if use_hinge:
            z = y_pred_or_z
            loss = np.mean(np.maximum(0, 1 - y_true * z))
        else:
            y_pred = np.clip(y_pred_or_z, eps, 1 - eps)
            loss = -np.mean(y_true * np.log(y_pred) + (1 - y_true) * np.log(1 - y_pred))
            
        if l2_lambda > 0 and weights is not None:
            loss += (l2_lambda / 2) * np.sum(weights**2)
        return loss

    def fit(self, X_train, y_train, X_val, y_val, epochs=100, lr=0.1, batch_size=32, 
            l2_lambda=0.0, use_momentum=False, momentum_beta=0.9, loss_type="bce"):
        """Обучение"""
        self._init_params(X_train.shape[1])
        self.train_loss_hist, self.val_loss_hist = [], []
        
        v_w, v_b = np.zeros_like(self.w), 0.0
        
        use_hinge = (loss_type == "hinge")
        y_train_h = np.where(y_train == 0, -1, 1) if use_hinge else y_train
        y_val_h = np.where(y_val == 0, -1, 1) if use_hinge else y_val

        n_samples = X_train.shape[0]
        for epoch in range(epochs):
            idx = np.random.permutation(n_samples)
            X_sh, y_sh = X_train[idx], y_train_h[idx]
            
            epoch_loss, n_batches = 0, 0
            for i in range(0, n_samples, batch_size):
                X_b, y_b = X_sh[i:i+batch_size], y_sh[i:i+batch_size]
                m = len(y_b)
                
                logits = X_b @ self.w + self.b
                y_pred = self.sigmoid(logits) if not use_hinge else logits
                
                if use_hinge:
                    mask = (y_b * logits < 1).astype(float)
                    dw = -(1/m) * (X_b.T @ (y_b * mask))
                    db = -(1/m) * np.sum(y_b * mask)
                else:
                    dw = (1/m) * (X_b.T @ (y_pred - y_b))
                    db = (1/m) * np.sum(y_pred - y_b)
                
                if l2_lambda > 0:
                    dw += (l2_lambda / m) * self.w
                
                if use_momentum:
                    v_w = momentum_beta * v_w + lr * dw
                    v_b = momentum_beta * v_b + lr * db
                    self.w -= v_w
                    self.b -= v_b
                else:
                    self.w -= lr * dw
                    self.b -= lr * db
                    
                epoch_loss += self.compute_loss(y_b, logits if use_hinge else y_pred, 
                                                l2_lambda=l2_lambda, use_hinge=use_hinge, weights=self.w)
                n_batches += 1
                
            self.train_loss_hist.append(epoch_loss / n_batches)
            val_logits = X_val @ self.w + self.b
            val_pred = self.sigmoid(val_logits) if not use_hinge else val_logits
            self.val_loss_hist.append(self.compute_loss(y_val_h, val_logits if use_hinge else val_pred, 
                                                        l2_lambda=l2_lambda, use_hinge=use_hinge, weights=self.w))
            
        return self

    def predict(self, X):
        return (self.forward(X) >= 0.5).astype(int)

    def predict_proba(self, X):
        return self.forward(X)