import React from "react";
import styles from "./Layout.scss";

const about = ABOUT || HOME || "/";

const Header = ({ tabs, onTabClick }) => {
  let children = tabs.map((tab, i) => {
    let css = styles.main_header_tab;
    if (tab.active) {
      css += " " + styles.active;
    }
    return (
      <span className={css} key={i} onClick={() => onTabClick(tab.id)}>
        <h2>{tab.id}</h2>
      </span>
    );
  });
  return (
    <div className={styles.main_header}>
      {children}
      <a href={about}>
        <span className={styles.main_header_tab} key="home">
          <h2>About</h2>
        </span>
      </a>
    </div>
  );
};

export default Header;
