import os
try:
    from PyQt6.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QLineEdit, QComboBox, 
                                 QPushButton, QHBoxLayout, QMessageBox, QFileDialog, QWidget, QCheckBox)
    from PyQt6.QtCore import Qt
except ModuleNotFoundError:
    from PyQt5.QtWidgets import (QDialog, QVBoxLayout, QFormLayout, QLineEdit, QComboBox, 
                                 QPushButton, QHBoxLayout, QMessageBox, QFileDialog, QWidget, QCheckBox)
    from PyQt5.QtCore import Qt

class SettingsDialog(QDialog):
    def __init__(self, config, parent=None):
        super().__init__(parent)
        self.config = config
        self.setWindowTitle("Connection Settings")
        self.setMinimumWidth(450)
        
        self.layout = QVBoxLayout(self)
        self.form_layout = QFormLayout()
        
        # Backend selection
        self.backend_cb = QComboBox()
        self.backend_cb.addItems(["local", "remote"])
        current_backend = self.config.get("Compute Backend", "local")
        self.backend_cb.setCurrentText(current_backend)
        self.form_layout.addRow("Compute Backend:", self.backend_cb)
        
        # URL
        self.url_le = QLineEdit(self.config.get("Remote Server URL", "https://localhost:8443"))
        self.form_layout.addRow("Remote Server URL:", self.url_le)
        
        # API Token
        self.token_le = QLineEdit(self.config.get("Remote API Token", ""))
        self.token_le.setEchoMode(QLineEdit.EchoMode.Password)
        self.form_layout.addRow("Remote API Token:", self.token_le)
        
        # Verify SSL
        self.verify_ssl_cb = QCheckBox("Verify SSL Certificates")
        self.verify_ssl_cb.setChecked(self.config.get("Remote Verify SSL") is not False)
        self.form_layout.addRow("", self.verify_ssl_cb)
        
        # Cert paths with browse buttons
        self.client_cert_le = self._add_file_picker("Remote Client Cert:", self.config.get("Remote Client Cert", ""))
        self.client_key_le = self._add_file_picker("Remote Client Key:", self.config.get("Remote Client Key", ""))
        self.ca_cert_le = self._add_file_picker("Remote CA Cert:", self.config.get("Remote CA Cert", ""))
        
        self.layout.addLayout(self.form_layout)
        
        # Buttons
        self.button_layout = QHBoxLayout()
        self.save_btn = QPushButton("Save && Restart Required")
        self.save_btn.clicked.connect(self.save_settings)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.clicked.connect(self.reject)
        
        self.button_layout.addWidget(self.save_btn)
        self.button_layout.addWidget(self.cancel_btn)
        self.layout.addLayout(self.button_layout)
        
        self._toggle_remote_fields()
        self.backend_cb.currentTextChanged.connect(self._toggle_remote_fields)

    def _add_file_picker(self, label, value):
        widget = QWidget()
        hlayout = QHBoxLayout(widget)
        hlayout.setContentsMargins(0, 0, 0, 0)
        
        le = QLineEdit(value)
        btn = QPushButton("...")
        btn.setMaximumWidth(30)
        
        def browse():
            file, _ = QFileDialog.getOpenFileName(self, f"Select {label.strip(':')}", os.path.dirname(le.text()) or ".", "All Files (*)")
            if file:
                le.setText(file)
                
        btn.clicked.connect(browse)
        
        hlayout.addWidget(le)
        hlayout.addWidget(btn)
        self.form_layout.addRow(label, widget)
        return le

    def _toggle_remote_fields(self):
        is_remote = self.backend_cb.currentText() == "remote"
        self.url_le.setEnabled(is_remote)
        self.token_le.setEnabled(is_remote)
        self.verify_ssl_cb.setEnabled(is_remote)
        self.client_cert_le.parent().setEnabled(is_remote)
        self.client_key_le.parent().setEnabled(is_remote)
        self.ca_cert_le.parent().setEnabled(is_remote)

    def save_settings(self):
        self.config.set("Compute Backend", self.backend_cb.currentText())
        self.config.set("Remote Server URL", self.url_le.text())
        self.config.set("Remote API Token", self.token_le.text())
        self.config.set("Remote Verify SSL", self.verify_ssl_cb.isChecked())
        self.config.set("Remote Client Cert", self.client_cert_le.text())
        self.config.set("Remote Client Key", self.client_key_le.text())
        self.config.set("Remote CA Cert", self.ca_cert_le.text())
        
        self.config.saveSettings()
        
        QMessageBox.information(self, "Settings Saved", "Connection settings have been saved. Please restart the application for them to take effect.")
        self.accept()
