def phase3(Sw,So,Sg,model="Stone's Model I",n=None):

    if self.Som is None:
        self._estimate_Som(Sg)

    if model=="Stone's Model I":
        self.kro,self.krw,self.krg = self._stones_model_I(Sw,So,Sg)
    elif model=="Aziz and Settari":
        self.kro,self.krw,self.krg = self._aziz_settari(Sw,So,Sg)
    elif model=="Stone's Model II":
        self.kro,self.krw,self.krg = self._stones_model_II(Sw,So,Sg)
    elif model=="Hustad-Holt Correlation":
        self.kro,self.krw,self.krg = self._hustad_holt(Sw,So,Sg,n)