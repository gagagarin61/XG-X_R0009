// ================================================
// X線画像 欠陥自動検出マクロ v1.0
// ================================================
// 検出対象:
//   [赤]    点欠損   - 孤立した異常ピクセル・クラスター
//   [シアン] 線欠損   - 水平・垂直方向の欠損ライン
//   [黄]    画像異常  - 局所的な輝度異常領域
// ------------------------------------------------
// 使い方:
//   1. X線画像を ImageJ で開く
//   2. このマクロを実行 (Plugins > Macros > Run...)
//   3. ダイアログでパラメータを調整して OK
//   4. 「欠陥マップ」画像と ROIマネージャーで結果を確認
// ================================================

macro "X線欠陥検出 [F2]" {

    if (nImages == 0) {
        exit("X線画像を先に開いてください。");
    }

    // ---- パラメータダイアログ ----
    Dialog.create("X線欠陥検出 設定");
    Dialog.addMessage("─── 検出感度 (σ値が大きいほど検出が少なくなる) ───");
    Dialog.addNumber("点欠損 閾値係数 σ:", 5.0);
    Dialog.addNumber("線欠損 閾値係数 σ:", 3.5);
    Dialog.addNumber("画像異常 閾値係数 σ:", 3.0);
    Dialog.addMessage("─── フィルタ設定 ───");
    Dialog.addNumber("点欠損 メディアンフィルタ半径 (px):", 2);
    Dialog.addNumber("画像異常 背景推定ブラー半径 (px):", 30);
    Dialog.addMessage("─── 検出項目の選択 ───");
    Dialog.addCheckbox("点欠損を検出する", true);
    Dialog.addCheckbox("線欠損を検出する", true);
    Dialog.addCheckbox("画像異常を検出する", true);
    Dialog.show();

    ptSigma  = Dialog.getNumber();
    lnSigma  = Dialog.getNumber();
    anSigma  = Dialog.getNumber();
    ptRadius = Dialog.getNumber();
    anBlur   = Dialog.getNumber();
    doPoint  = Dialog.getCheckbox();
    doLine   = Dialog.getCheckbox();
    doAnom   = Dialog.getCheckbox();

    // ---- 元画像情報 ----
    origID    = getImageID();
    origTitle = getTitle();
    W = getWidth();
    H = getHeight();

    setBatchMode(true);

    // 32bit グレースケール作業コピー
    selectImage(origID);
    run("Duplicate...", "title=_xray_work");
    gID = getImageID();
    if (bitDepth() != 32) run("32-bit");

    roiManager("reset");
    ptCount = 0;
    lnCount = 0;
    anCount = 0;

    // ================================================
    // A. 点欠損検出 (Point Defect Detection)
    //    原理: メディアンフィルタとの差分で孤立した異常ピクセルを検出
    // ================================================
    if (doPoint) {
        showStatus("点欠損を検出中...");
        selectImage(gID);

        // メディアンフィルタで周辺画素の期待値を推定
        run("Duplicate...", "title=_pt_med");
        ptMedID = getImageID();
        run("Median...", "radius=" + ptRadius);

        // 差分 = 原画像 - メディアン → 突出した画素が大きな値になる
        selectImage(gID);
        imageCalculator("Subtract create 32-bit", gID, ptMedID);
        ptDiffID = getImageID();

        // 絶対値化して明・暗両方向の欠損を検出
        run("Abs");
        getStatistics(area, ptMean, ptMin, ptMax, ptStd);
        ptThresh = ptMean + ptSigma * ptStd;

        // 閾値でマスク化
        setThreshold(ptThresh, 1e15);
        setOption("BlackBackground", true);
        run("Convert to Mask");

        // 連結成分解析で孤立した小領域を検出 (size=1-500 px)
        ptROIStart = roiManager("count");
        run("Analyze Particles...", "size=1-500 circularity=0.0-1.0 display clear add");
        ptCount = nResults;

        // ROIを赤でマーク
        for (i = 0; i < ptCount; i++) {
            roiManager("select", ptROIStart + i);
            Roi.setStrokeColor("red");
            Roi.setStrokeWidth(1);
            roiManager("update");
        }

        selectImage(ptMedID);  close();
        selectImage(ptDiffID); close();
    }

    // ================================================
    // B. 線欠損検出 (Line Defect Detection)
    //    原理: 行・列ごとの平均輝度を計算し、統計的外れ値の行・列を検出
    // ================================================
    if (doLine) {
        showStatus("線欠損を検出中...");
        selectImage(gID);

        rowMeans = newArray(H);
        colMeans = newArray(W);

        // 各行の平均輝度を取得
        for (y = 0; y < H; y++) {
            showProgress(y, H + W);
            makeRectangle(0, y, W, 1);
            getStatistics(a, m);
            rowMeans[y] = m;
        }

        // 各列の平均輝度を取得
        for (x = 0; x < W; x++) {
            showProgress(H + x, H + W);
            makeRectangle(x, 0, 1, H);
            getStatistics(a, m);
            colMeans[x] = m;
        }

        run("Select None");

        // 行平均の平均と標準偏差を計算
        rSum = 0; rSum2 = 0;
        for (y = 0; y < H; y++) {
            rSum  += rowMeans[y];
            rSum2 += rowMeans[y] * rowMeans[y];
        }
        rMean = rSum / H;
        rStd  = sqrt(rSum2 / H - rMean * rMean);

        // 列平均の平均と標準偏差を計算
        cSum = 0; cSum2 = 0;
        for (x = 0; x < W; x++) {
            cSum  += colMeans[x];
            cSum2 += colMeans[x] * colMeans[x];
        }
        cMean = cSum / W;
        cStd  = sqrt(cSum2 / W - cMean * cMean);

        // 異常な行 (水平線欠損) を検出
        for (y = 0; y < H; y++) {
            if (abs(rowMeans[y] - rMean) > lnSigma * rStd) {
                makeRectangle(0, y, W, 1);
                roiManager("add");
                roiManager("select", roiManager("count") - 1);
                roiManager("rename", "H-Line_y=" + y);
                Roi.setStrokeColor("cyan");
                Roi.setStrokeWidth(1);
                roiManager("update");
                lnCount++;
            }
        }

        // 異常な列 (垂直線欠損) を検出
        for (x = 0; x < W; x++) {
            if (abs(colMeans[x] - cMean) > lnSigma * cStd) {
                makeRectangle(x, 0, 1, H);
                roiManager("add");
                roiManager("select", roiManager("count") - 1);
                roiManager("rename", "V-Line_x=" + x);
                Roi.setStrokeColor("cyan");
                Roi.setStrokeWidth(1);
                roiManager("update");
                lnCount++;
            }
        }

        run("Select None");
    }

    // ================================================
    // C. 画像異常検出 (Image Anomaly Detection)
    //    原理: 大きいガウシアンブラーで背景推定し、差分で局所異常を検出
    // ================================================
    if (doAnom) {
        showStatus("画像異常を検出中...");
        selectImage(gID);

        // 強いブラーで低周波背景を推定
        run("Duplicate...", "title=_anom_bg");
        anBGID = getImageID();
        run("Gaussian Blur...", "sigma=" + anBlur);

        // 差分 = 原画像 - 背景 → 局所的な輝度変動が残る
        selectImage(gID);
        imageCalculator("Subtract create 32-bit", gID, anBGID);
        anDiffID = getImageID();

        run("Abs");
        getStatistics(area, anMean, anMin, anMax, anStd);
        anThresh = anMean + anSigma * anStd;

        setThreshold(anThresh, 1e15);
        setOption("BlackBackground", true);
        run("Convert to Mask");

        // 中程度以上の面積を持つ領域を検出 (size=100px以上)
        anROIStart = roiManager("count");
        run("Analyze Particles...", "size=100-Infinity display clear add");
        anCount = nResults;

        // ROIを黄色でマーク
        for (i = 0; i < anCount; i++) {
            roiManager("select", anROIStart + i);
            Roi.setStrokeColor("yellow");
            Roi.setStrokeWidth(2);
            roiManager("update");
        }

        selectImage(anBGID);   close();
        selectImage(anDiffID); close();
    }

    // ================================================
    // 結果画像の作成 (ROIオーバーレイ付き)
    // ================================================
    selectImage(origID);
    run("Duplicate...", "title=[欠陥マップ: " + origTitle + "]");
    resID = getImageID();
    if (bitDepth() != 24) run("RGB Color");

    selectImage(gID); close();

    setBatchMode("exit and display");

    // 結果画像にROIオーバーレイを表示
    selectImage(resID);
    if (roiManager("count") > 0) {
        roiManager("show all without labels");
        // 画像にROIを焼き付け
        run("Flatten");
        flatID = getImageID();
        rename("欠陥マップ: " + origTitle);
        selectImage(resID); close();
    }

    // ================================================
    // ROIマネージャーとログに結果を出力
    // ================================================
    total = ptCount + lnCount + anCount;

    // ログウィンドウに詳細を出力
    print("\\Clear");
    print("================================================");
    print(" X線欠陥検出結果");
    print("================================================");
    print(" 対象画像 : " + origTitle);
    print(" サイズ   : " + W + " x " + H + " px");
    print("------------------------------------------------");
    if (doPoint) print(" [赤]    点欠損   : " + ptCount + " 箇所");
    if (doLine)  print(" [シアン] 線欠損   : " + lnCount + " ライン");
    if (doAnom)  print(" [黄]    画像異常  : " + anCount + " 領域");
    print("------------------------------------------------");
    print(" 合計 : " + total + " 件");
    print("================================================");
    print(" ※ ROIマネージャーで個別選択・計測が可能です");
    if (doPoint) print(" ※ 点欠損の詳細は Results テーブルを参照");

    // サマリーダイアログ
    msg  = "=== X線欠陥検出完了 ===\n\n";
    msg += "対象: " + origTitle + "\n";
    msg += "サイズ: " + W + " × " + H + " px\n\n";
    if (doPoint) msg += "■ 点欠損   [赤]    : " + ptCount + " 箇所\n";
    if (doLine)  msg += "■ 線欠損   [シアン]: " + lnCount + " ライン\n";
    if (doAnom)  msg += "■ 画像異常 [黄]    : " + anCount + " 領域\n";
    msg += "\n合計: " + total + " 件検出\n\n";
    msg += "・欠陥マップ画像に結果が表示されます\n";
    msg += "・ROIマネージャーで個別選択・計測が可能です\n";
    msg += "・詳細はログウィンドウを参照してください";

    showMessage("検出結果", msg);
}


// ================================================
// サブマクロ: 点欠損のみ再検出
// ================================================
macro "点欠損のみ検出 [F3]" {
    if (nImages == 0) exit("画像を開いてください。");

    ptSigma  = getNumber("点欠損 閾値係数 σ:", 5.0);
    ptRadius = getNumber("メディアンフィルタ半径 (px):", 2);

    origID = getImageID();
    setBatchMode(true);

    run("Duplicate...", "title=_pt_tmp");
    tmpID = getImageID();
    run("32-bit");

    run("Duplicate...", "title=_pt_med");
    medID = getImageID();
    run("Median...", "radius=" + ptRadius);

    imageCalculator("Subtract create 32-bit", tmpID, medID);
    diffID = getImageID();
    run("Abs");

    getStatistics(a, m, mn, mx, s);
    setThreshold(m + ptSigma * s, 1e15);
    setOption("BlackBackground", true);
    run("Convert to Mask");

    roiManager("reset");
    run("Analyze Particles...", "size=1-500 display clear add");
    count = nResults;

    for (i = 0; i < count; i++) {
        roiManager("select", i);
        Roi.setStrokeColor("red");
        roiManager("update");
    }

    selectImage(tmpID);  close();
    selectImage(medID);  close();
    selectImage(diffID); close();

    setBatchMode("exit and display");
    selectImage(origID);
    roiManager("show all without labels");

    showMessage("点欠損検出", count + " 箇所の点欠損を検出しました。");
}


// ================================================
// サブマクロ: 線欠損のみ再検出
// ================================================
macro "線欠損のみ検出 [F4]" {
    if (nImages == 0) exit("画像を開いてください。");

    lnSigma = getNumber("線欠損 閾値係数 σ:", 3.5);

    origID = getImageID();
    W = getWidth();
    H = getHeight();

    setBatchMode(true);

    run("Duplicate...", "title=_ln_tmp");
    tmpID = getImageID();
    run("32-bit");

    rowMeans = newArray(H);
    colMeans = newArray(W);

    for (y = 0; y < H; y++) {
        makeRectangle(0, y, W, 1);
        getStatistics(a, m);
        rowMeans[y] = m;
    }
    for (x = 0; x < W; x++) {
        makeRectangle(x, 0, 1, H);
        getStatistics(a, m);
        colMeans[x] = m;
    }
    run("Select None");

    rSum = 0; rSum2 = 0;
    for (y = 0; y < H; y++) { rSum += rowMeans[y]; rSum2 += rowMeans[y]*rowMeans[y]; }
    rMean = rSum/H; rStd = sqrt(rSum2/H - rMean*rMean);

    cSum = 0; cSum2 = 0;
    for (x = 0; x < W; x++) { cSum += colMeans[x]; cSum2 += colMeans[x]*colMeans[x]; }
    cMean = cSum/W; cStd = sqrt(cSum2/W - cMean*cMean);

    roiManager("reset");
    lnCount = 0;

    for (y = 0; y < H; y++) {
        if (abs(rowMeans[y] - rMean) > lnSigma * rStd) {
            makeRectangle(0, y, W, 1);
            roiManager("add");
            roiManager("select", roiManager("count") - 1);
            roiManager("rename", "H-Line_y=" + y);
            Roi.setStrokeColor("cyan");
            roiManager("update");
            lnCount++;
        }
    }
    for (x = 0; x < W; x++) {
        if (abs(colMeans[x] - cMean) > lnSigma * cStd) {
            makeRectangle(x, 0, 1, H);
            roiManager("add");
            roiManager("select", roiManager("count") - 1);
            roiManager("rename", "V-Line_x=" + x);
            Roi.setStrokeColor("cyan");
            roiManager("update");
            lnCount++;
        }
    }
    run("Select None");

    selectImage(tmpID); close();

    setBatchMode("exit and display");
    selectImage(origID);
    roiManager("show all without labels");

    showMessage("線欠損検出", lnCount + " ラインの線欠損を検出しました。\n(水平: 行の異常, 垂直: 列の異常)");
}
